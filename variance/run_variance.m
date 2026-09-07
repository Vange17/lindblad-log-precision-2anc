function run_variance(num_NMt, seed, output_stem)
% Sequential variance pilot using the pinned upstream sampling loop.
% Only scientific correction: initialize xlab_prod before every kA branch.
% Equivalent implementation changes: cached UH; Kraus partial trace;
% local inverse-CDF categorical sampler; no symbolic preprocessing or shots.
validateattributes(num_NMt,{'numeric'},{'scalar','integer','>=',2,'<=',10000});
validateattributes(seed,{'numeric'},{'scalar','integer','>=',0,'<=',2^32-1});
assert(ischar(output_stem) && isrow(output_stem),'Provide an output file stem.');
assert(~isfile([output_stem '.mat']) && ~isfile([output_stem '.csv']),...
 'Output already exists; use a new output stem.');
output_dir=fileparts(output_stem);
assert(isempty(output_dir) || isfolder(output_dir),'Create the output directory first.');
variance_dir=fileparts(mfilename('fullpath'));
mat_path=fullfile(fileparts(variance_dir),...
 'miti_Ising_OQS4_sp5n500gam1p5rho1O1_NMnoshots.mat');
load(mat_path,'Ham','Jl','rho0','rho11','num_spin','dt','t_vec',...
 'A','cla1','PL','p_cla1','cHa0','PH','xl','mundt1','mumdt1',...
 'sumeqk','p_ss','mean_overlap_rho0J0','overlap_rhoexact_mean',...
 'mean_overlap_rhoNM');
rng(seed,'twister');
transA=@(x) dec2base(x,3)-'0';
UH=expm(-1i*Ham*dt);
d=2^num_spin;
K0=Jl(1:d,1:d); K1=Jl(d+1:2*d,1:d);
baseline=zeros(size(t_vec)); rho=rho0;
for j=1:numel(t_vec)
 if j>1, rho=UH*rho*UH'; rho=K0*rho*K0'+K1*rho*K1'; end
 baseline(j)=real(trace(rho11*rho));
end
baseline_max_difference=max(abs(baseline-mean_overlap_rho0J0));
assert(baseline_max_difference<1e-11);
mu_step=sum(mundt1)*sum(mumdt1);
overlap_rhoNM=zeros(num_NMt,numel(t_vec));
timer=tic;
for qNMt=1:num_NMt
    count_tt=0;
    overlap_rhoNM_tem=zeros(1,length(t_vec)); %预分配内存

    for tt=t_vec  %tt必须从第一个dt开始，因为后面作用的Jl里面包含了一个dt
        count_tt=count_tt+1;
        if count_tt==1
            rho0_Hdt=rho0;
            rho_NM=rho0;
            overlap_rhoNM_tem(count_tt)=trace(rho11*rho0);
            continue;

        else
            rho0_Hdt=UH*rho_NM*UH';
        end

        rho0_J0=K0*rho0_Hdt*K0'+K1*rho0_Hdt*K1'; %注意：这里的rho0_J0其实是t时刻时经过随机mitigation之后的rho_NM，此处rho_NM未必是合法量子态，未必厄米。

        p_kmundt1=mundt1/sum(mundt1);
        k_mundt1=randsrc(1,1,[1:length(mundt1);p_kmundt1]); %生成k阶次
        Pa=eye(2^num_spin); Pb=eye(2^num_spin);

        kA=k_mundt1-1;
        xlab_prod=1; % Local correction: reset at every time step.
        if kA==0
            Pa;Pb;
        else
            %---------------------------%↓生成A{kA}的抽样及相关编号------------------------------
            pA=abs(A{kA})/sum(abs(A{kA}));
            length_A = 1:length(A{kA});
            pos_A=randsrc(1,1,[length_A;pA]);
            samples_A=A{kA}(pos_A);

            pos_A3=transA(pos_A-1); %转换成3进制数字编号.与trans1tombit_add1函数共用时，千万别忘记有可能需要-1处理！
            pos_A3add1=trans1tombit_add1(pos_A3,kA,3);%转成{1,2,3}构成的3进制形式
            %---------------------------%↑生成A{kA}的抽样及相关编号------------------------------

            %----------%↓生成chi_ab^{1}chi_ab^{2}...chi_ab^{k}连乘以及P1P2 ρ P2P3...的连乘--------
            xl_ab_tem=1;
% %             rho_tem=rho0_J0; %重置变量
            for k=kA:-1:1%k表示采样到A^{k}阶  %注意：当for循环降序循环时，中间的值-1必不可少！！！
                class_ADi=pos_A3add1(k); %class_ADi表抽样到的ji编号，ji∈{1,2,3}
                a_xlab=randsrc(1,1,[1:length(cla1);p_cla1]);b_xlab=randsrc(1,1,[1:length(cla1);p_cla1]);
                cla1_randa=cla1(a_xlab); cla1_randb=cla1(b_xlab);
                xl_ab_tem=xl_ab_tem*cla1_randa*cla1_randb';
                % %                             xl_ab_tem=exp( 1i*angle(cla1_randa*cla1_randb') );
                if class_ADi==1  %P1P2 ρ P2P3...的连乘的构型判断
                    xl_ab_tem;
                    Pa=PL{a_xlab}*Pa;Pb=PL{b_xlab}*Pb;
                elseif class_ADi==2
                    xl_ab_tem;%注：此处不是xl_ab_tem=sign(-1/2)*xl_ab_tem，因为正负号信息已经包含在后续对应的A_j1j2...jkA里了！
                    Pa=PL{b_xlab}'*PL{a_xlab}*Pa;Pb=eye(2^num_spin)*Pb;
                else             %==3
                    xl_ab_tem;
                    Pa=eye(2^num_spin)*Pa;Pb=(PL{b_xlab}'*PL{a_xlab})'*Pb;
                end
            end
            xlab_prod=xlab_prod*A{kA}( transkAto10(pos_A3add1) )*xl_ab_tem;         %chi_ab^{1}chi_ab^{2}...chi_ab^{k}的连乘结果
            %----------%↑生成chi_ab^{1}chi_ab^{2}...chi_ab^{k}连乘以及P1P2 ρ P2P3...的连乘--------
            %----------------------------------------------------------------------------

        end


        p_kmumdt1=mumdt1/sum(mumdt1);
        pos_mumdt1=randsrc(1,1,[1:length(mumdt1);p_kmumdt1]);

        kM=pos_mumdt1-1; %一共C_(kM+2)^2个。kM表dt^kM项

        if kM==0 %[s s0 s1...]采样程序，储存结果为samples_ss。
            ss=sumeqk{kM+1}; %ss=[0 0 0]
            pos_ss=1;
            samples_ss=ss;
        else
            ss=sumeqk{kM+1};    %生成[s s0 s1...sm]序列,要求sum(s)=k。此例m=1.%sumeqk事先已生成好，此处仅查询以提高效率。%一共C_(kM+2)^2个
            pos_ss=randsrc(1,1,[1:size(ss, 1);p_ss{kM+1}]);
            samples_ss=ss(pos_ss,:);
        end
        samples_ss;
        %----------s1,s2,...sm段信息↓-----------------
        if samples_ss(3)==0
            xlab_prod=1*xlab_prod;
        else %samples_ss(3)>=1
            xlab_prod=1*xlab_prod;
            for qsm=samples_ss(3):-1:1 %qsm表sm段的编号
                %产生xl_ab编号及PρP...操作
                class_xlab=randsrc(1,1,[1:3;1/2 1/4 1/4]); %从xl_ab的3种构型中选择1种
                a_xlab=randsrc(1,1,[1:length(cla1);p_cla1]);b_xlab=randsrc(1,1,[1:length(cla1);p_cla1]);
                cla1_randa=cla1(a_xlab); cla1_randb=cla1(b_xlab);
                xlab_prod= xlab_prod*cla1_randa*cla1_randb'; %从具体的10进制编号中获取xlab信息,连乘
                if class_xlab==1
                    xlab_prod;
                    Pa=PL{a_xlab}*Pa;Pb=PL{b_xlab}*Pb;
                elseif class_xlab==2
                    xlab_prod=sign(-1/2)*xlab_prod;
                    Pa=PL{b_xlab}'*PL{a_xlab}*Pa;Pb=eye(2^num_spin)*Pb;
                else %class_xlab==3
                    xlab_prod=sign(-1/2)*xlab_prod;
                    Pa=eye(2^num_spin)*Pa;Pb=(PL{b_xlab}'*PL{a_xlab})'*Pb;
                end
            end
        end
        %----------s1,s2,...sm段信息↑-----------------

        %----------s0段信息↓-----------------
        if samples_ss(2)==0 %各变量保持不变
            xlab_prod;
            Pa;Pb;
% %             rho_Mtem;
        else %samples_ss(2)>=1
            for qH=samples_ss(2):-1:1 %qH表H段的编号;%产生xl_ab编号及PρP...操作
                p_H=abs(cHa0)/sum(abs(cHa0));
                pos_H=randsrc(1,1,[1:length(cHa0);p_H]);
                class_H=randsrc(1,1,[1:2;1/2 1/2]);%超算子H共2中构型，分别是-iHρ与(-iHρ)',两种构型出现概率相同！
                if class_H==1
                    xlab_prod=-1i*cHa0(pos_H)*xlab_prod;
                    Pa=PH{pos_H}*Pa;Pb=eye(2^num_spin)*Pb;
                else %class_H==2
                    xlab_prod=( -1i*cHa0(pos_H) )'*xlab_prod;
                    Pa=eye(2^num_spin)*Pa;Pb=PH{pos_H}*Pb;
                end
            end
        end
        %----------s0段信息↑-----------------

        %----------s段信息↓-----------------
        if samples_ss(1)==0
            xlab_prod;
            Pa;Pb;
        else %samples_ss(1)>=1
            for qHsm=samples_ss(1):-1:1 %qH表H段的编号;%产生xl_ab编号及PρP...操作
                Prji_s=xl/sum(xl);
                class_Hsm=randsrc(1,1,[1:length(xl);Prji_s]); %产生s位ji信息，ji=0,1,...m。
                if class_Hsm==1     %ji_s表s段的编号对应的信息值
                    p_H=abs(cHa0)/sum(abs(cHa0));
                    pos_H=randsrc(1,1,[1:length(cHa0);p_H]);
                    class_H=randsrc(1,1,[1:2;1/2 1/2]);%超算子H共2中构型，分别是-iHρ与(-iHρ)',两种构型出现概率相同！
                    if class_H==1 %此编号装填的是xl0部分信息，即H超算子部分信息
                        xlab_prod=-1i*cHa0(pos_H)*xlab_prod;
                        Pa=PH{pos_H}*Pa;Pb=eye(2^num_spin)*Pb;
                    else %class_H==2
                        xlab_prod=(-1i*cHa0(pos_H))'*xlab_prod;
                        Pa=eye(2^num_spin)*Pa;Pb=PH{pos_H}*Pb;
                    end
                else %class_Hsm>=2。%此编号装填的是xl(2:end)部分信息，即L1...Lm耗散超算子部分信息
                    %产生xl_ab编号及PρP...操作
                    class_xlab=randsrc(1,1,[1:3;1/2 1/4 1/4]); %从xl_ab的3种构型中选择1种
                    a_xlab=randsrc(1,1,[1:length(cla1);p_cla1]);b_xlab=randsrc(1,1,[1:length(cla1);p_cla1]);
                    cla1_randa=cla1(a_xlab); cla1_randb=cla1(b_xlab);
                    xlab_prod= xlab_prod*cla1_randa*cla1_randb'; %从具体的10进制编号中获取xlab信息,连乘
                    if class_xlab==1
                        xlab_prod;
                        Pa=PL{a_xlab}*Pa;Pb=PL{b_xlab}*Pb;
                    elseif class_xlab==2
                        xlab_prod=sign(-1/2)*xlab_prod;
                        Pa=PL{b_xlab}'*PL{a_xlab}*Pa;Pb=eye(2^num_spin)*Pb;
                    else %class_xlab==3
                        xlab_prod=sign(-1/2)*xlab_prod;
                        Pa=eye(2^num_spin)*Pa;Pb=(PL{b_xlab}'*PL{a_xlab})'*Pb;
                    end
                end
            end
        end
        %----------s段信息↑-----------------
        xlab_prod=xlab_prod*(-1)^( kM-samples_ss(1) );%之所以有该行，是因为M前有系数(-1)^(k-s)。kM表dt^kM项
        phi=angle(xlab_prod);
        rho_M_tem=exp(1i*phi)*Pa*rho0_J0*Pb'; %做shots时不能乘以系数sum(mumdt1)，因为有可能使overlap>1，使shots失败
        rho_M=(rho_M_tem+rho_M_tem')/2;

        rho_NM=rho_M;

        overlap_rhoNM_tem(count_tt)=trace(rho11*rho_NM); %加abs是为了防止出现特别接近于0的负数出现！
    end
    overlap_rhoNM(qNMt,:)=overlap_rhoNM_tem;%overlap_rhoNM实际上是观测量O的均值！只是此时恰好O=|1><1|，均值刚好是overlap.
end


elapsed_seconds=toc(timer);
values=real(overlap_rhoNM).*(mu_step.^(0:numel(t_vec)-1));
assert(max(abs(imag(overlap_rhoNM)),[],'all')<1e-10);
trajectory_mean=mean(values,1);
trajectory_variance=var(values,0,1);
trajectory_std=sqrt(trajectory_variance);
standard_error_this_run=sqrt(trajectory_variance/num_NMt);
standard_error_N5million=sqrt(trajectory_variance/5e6);
exact_mean=real(overlap_rhoexact_mean);
archive_mean=mean_overlap_rhoNM;
software_version=version;
save([output_stem '.mat'],'num_NMt','seed','t_vec','trajectory_mean',...
 'trajectory_variance','trajectory_std','standard_error_this_run','standard_error_N5million',...
 'exact_mean','archive_mean','mu_step','elapsed_seconds','baseline_max_difference',...
 'software_version','values');
output=table(t_vec(:),trajectory_mean(:),trajectory_variance(:),trajectory_std(:),...
 standard_error_this_run(:),standard_error_N5million(:),exact_mean(:),archive_mean(:),...
 'VariableNames',{'time','mean','variance','trajectory_std','se_this_run','se_N5million','exact_mean','archive_mean'});
writetable(output,[output_stem '.csv']);
fprintf('Completed N=%d, seed=%d in %.3f seconds; baseline difference %.3g\n',...
 num_NMt,seed,elapsed_seconds,baseline_max_difference);
disp(output([1 26 51 76],:));
end

function value=randsrc(rows,cols,distribution)
assert(rows==1 && cols==1);
probabilities=distribution(2,:);
assert(all(probabilities>=0) && abs(sum(probabilities)-1)<1e-10);
cdf=cumsum(probabilities); cdf(end)=1;
value=distribution(1,find(rand<cdf,1));
end

function result=trans1tombit_add1(mat,mbit,nary)
assert(isvector(mat) && size(mat,1)==1 && numel(mat)<=mbit && nary==3);
result=[ones(1,mbit-numel(mat)),mat+1];
end

function result=transkAto10(mat)
result=1+sum((mat-1).*3.^(numel(mat)-1:-1:0));
end

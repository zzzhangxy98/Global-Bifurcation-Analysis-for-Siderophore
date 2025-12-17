spall=cell(1,32);
parfor i=1:32
    R_supply = 1;
    a_11=0.3; a_22=0.4; a_33=0.5;
    alpha_matrix = [a_11,0,0; 0,a_22,0; 0,0,a_33];
    v123 = 0.64+0.005*i;
    v321 = 0;
    v_matrix = [1-v123, v123, v321;
                v321, 1-v123, v123;
                v123, v321, 1-v123];
    opt=struct;
    opt.d=0.2;
    opt.dt1=0.002;
    opt.dt2=0.002;
    opt.max=200000;
    opt.m=1;
    opt.l=1e-4;
    opt.eps=1e-8;
    opt.s=1;
    opt.n=7;
    opt.alpha_sid_pro = alpha_matrix;
    opt.v_sid_rec = v_matrix;
    opt.e = 10;
    opt.u = 1;
    opt.migr = 0;
    opt.gamma = 1;
    opt.d0 = 0.1;
    opt.R_sup = R_supply;
    F_func=@ODE_siderophore_iron_partition;
    
    
        
    [sp]=makesl4TIEZAIT(F_func,num2str(i),opt)
    spall{i}=sp;
end
save(['spall','_ODE_siderophore_iron_partition','.mat'],"spall")
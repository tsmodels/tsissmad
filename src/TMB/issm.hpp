/// @file issm.hpp
#ifndef issm_hpp
#define issm_hpp

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

template<class Type>
Type issm(objective_function<Type>* obj) {
    DATA_VECTOR(V);
    DATA_MATRIX(X);
    DATA_VECTOR(good);
    DATA_VECTOR(y);
    DATA_VECTOR(allpars);
    DATA_IVECTOR(findex);
    DATA_IVECTOR(fpindex);
    DATA_IVECTOR(ppindex);
    DATA_IVECTOR(fshape);
    DATA_IVECTOR(modeli);
    PARAMETER_VECTOR(pars);
    PARAMETER(lambda);
    PARAMETER_VECTOR(kappa);
    // slot the parameters into the vectors
    for(int i=0;i<ppindex.size();i++){
        allpars(ppindex(i)) = pars(i);
    }
    for(int i = 0;i<findex.size();i++){
        V(findex(i)) = allpars(fpindex(i));
    }
    vector<Type> F0tmp = V.segment(fshape(0),fshape(1));
    matrix<Type> F0 = asMatrix(F0tmp, modeli(0),modeli(0));
    vector<Type> F1tmp = V.segment(fshape(2),fshape(3));
    matrix<Type> F1 = asMatrix(F1tmp, modeli(0),modeli(0));
    vector<Type> F2tmp = V.segment(fshape(4),fshape(5));
    matrix<Type> F2 = asMatrix(F2tmp, modeli(0),modeli(0));
    matrix<Type> F = (F0.array() * F1.array() * F2.array()).matrix();
    vector<Type> xreg = X * kappa;
    vector<Type> W = V.segment(fshape(6),fshape(7));
    vector<Type> G = V.segment(fshape(8),fshape(9));
    int timesteps = y.size();
    matrix<Type> xaux(timesteps, modeli(0));
    matrix<Type> waux(timesteps, modeli(0));
    matrix<Type> states(timesteps + 1, modeli(0));
    vector<Type> eaux(timesteps);
    vector<Type> error(timesteps + 1);
    vector<Type> yaux(timesteps);
    vector<Type> yhat(timesteps + 1);
    xaux.setZero();
    waux.setZero();
    eaux.setZero();
    yaux.setZero();
    error.setZero();
    yhat.setZero();
    states.setZero();
    vector<Type> ytrans(y.size());
    for(int i = 0;i<timesteps;i++){
        ytrans(i) = CppAD::CondExpLe(lambda, Type(1.0e-11), log(y(i)), (pow(y(i), lambda) - Type(1.0))/lambda);
    }
    yaux(0) = xreg(0);
    waux.row(0) = W;
    eaux(0) = ytrans(0) - yaux(0);
    xaux.row(0) = G * eaux(0);
    matrix<Type> M = asMatrix(W,modeli(0),1).transpose();
    matrix<Type> B = asMatrix(G,modeli(0),1) * M;
    matrix<Type> D = F.array() -  B.array();
    vector<Type> tmp = xaux.row(0);
    vector<Type> gtmp = F * xaux.row(0).transpose();
    for (int i = 1; i < timesteps; i++) {
        tmp = xaux.row(i-1);
        yaux(i) = (tmp.array() * W.array()).sum() + xreg(i);
        if (good(i) > 0.5) {
            eaux(i) = ytrans(i) - yaux(i);
        } else {
            eaux(i) = 0.0;
        }
        gtmp = F * xaux.row(i-1).transpose();
        xaux.row(i) = (gtmp + G * eaux(i));
        waux.row(i) = waux.row(i-1) * D;
    }
    vector<Type> ytranst(timesteps + 1);
    vector<Type> goodt(timesteps + 1);
    vector<Type> xregt(timesteps + 1);
    
    matrix<Type> A = waux.leftCols(modeli(1));
    matrix<Type> init_states = A.householderQr().solve(eaux.matrix());
    states.topLeftCorner(1,modeli(1)) = init_states.transpose();
    Type sum_error_squared = 0.0;
    Type sum_log_y = 0.0;
    vector<Type> tmp2(modeli(0));
    vector<Type> gtmp2 = F * states.row(0).transpose();
    for(int i = 1; i <= timesteps; i++) {
        tmp2 = states.row(i-1);
        yhat(i) = (tmp2.array() * W.array()).sum() + xreg(i - 1);
        if(good(i - 1) > 0.5) {
            error(i) = ytrans(i - 1) - yhat(i);
            sum_error_squared+= error(i) * error(i);
            sum_log_y+=log(y(i-1));
        } else {
            error(i) = 0.0;
        }
        gtmp2 = F * states.row(i-1).transpose();
        states.row(i) = (gtmp2 + G * error(i));
    }
    Type good_timesteps = good.sum();
    REPORT(D);
    REPORT(states);
    Type loglik = (good_timesteps) * log(sum_error_squared) - Type(2.0) * (lambda - Type(1.0)) * sum_log_y;
    return(loglik);
}

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this

#endif
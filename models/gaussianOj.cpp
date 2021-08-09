#include <RcppArmadillo.h>
using namespace arma;
using namespace std;
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]

class gaussianObj{
public:
    
    // y_ : list of observed values
    // bMatList: list of 
    gaussianObj(mat Ymat_, mat Zmat_){
        Ymat = Ymat_;
        Zmat = Zmat_;
        nTotal = Ymat.n_rows;
        mTotal = Ymat.n_cols;
        YtZGlobal = Ymat.t()*Zmat/nTotal;
        ZtZGlobal = Zmat.t() * Zmat/nTotal;
        lambda = 0;
        
    }
    // we might have further speed improvment with refinement of bMatLarge?
    
    
    void set_penaltyMatrix(mat Omega_){
        Omega = Omega_;
    }
    
    void set_tuningParameter(double lambda_){
        lambda = lambda_;
    }
    
    
    // objective function for manifold
    // C = U * W * V^T. U\in R^{K\times r}. W\in R^{r\times r}
    double objF(List UWVt){
        mat U, W, V, Chat;
        mat yHat, yDiff;
        double loss;
        
        U = as<mat>(UWVt(0));
        W = as<mat>(UWVt(1));
        V = as<mat>(UWVt(2));
        Chat = U * W * V.t();
        yHat = Zmat * Chat;
        yDiff = Ymat - yHat;
        loss = accu(yDiff % yDiff)/nTotal;
        loss += lambda * accu(Chat % (Omega* Chat));
        return loss;
    }
    
    // gradient function for manifold
    // Input: List(U, W, V)
    // Output: List(gradU, gradW, gradV)
    List gradF(List UWVt){
        mat U, W, V, Chat;
        mat yHat, yDiff;
        mat gradU, gradW, gradV;

        U = as<mat>(UWVt(0));
        W = as<mat>(UWVt(1));
        V = as<mat>(UWVt(2));
        Chat = U * W * V.t();
        yHat = Zmat * Chat;
        yDiff = Ymat - yHat;
        
        mat tmp;
        tmp = ZtZGlobal+lambda*Omega;
        
        gradV = -2*YtZGlobal * U * W;
        gradU = -2 *YtZGlobal.t()*V*W + 2*tmp*U*W*W;
        
        tmp = U.t() * tmp * U * W;
        gradW = -2*U.t()*YtZGlobal.t()*V+tmp+tmp.t();
        
        List gradL = List::create(gradU, gradW, gradV);
        return gradL;
    }
    
    
private:
    // total obs points
    size_t nTotal;
    // total number of tasks
    size_t mTotal;
    
    //tuning parameters on rank and smoothness respectively
    double lambda;
    mat Ymat, Zmat;
    mat  Omega;
    mat YtZGlobal, ZtZGlobal;
};



RCPP_MODULE(gaussianObj_MODUDLE){
    class_<gaussianObj>("gaussianObj")
    .constructor<mat, mat>()
    .method("set_penaltyMatrix", &gaussianObj::set_penaltyMatrix)
    .method("set_tuningParameter", &gaussianObj::set_tuningParameter)
    .method("objF", &gaussianObj::objF)
    .method("gradF", &gaussianObj::gradF)
    ;

}

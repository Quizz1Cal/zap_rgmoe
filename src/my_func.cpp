#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
List my_func() {

    CharacterVector x = CharacterVector::create( "foo", "bar" )  ;
    NumericVector y   = NumericVector::create( 0.0, 1.0 ) ;
    List z            = List::create( x, y ) ;

    return z ;
}

// [[Rcpp::export]]
LogicalVector all(LogicalVector b) {
    int n = b.size();
    for (int i=0; i < n; i++) {
        if (!b[i]) {
            return false;
        }
    }
    return true;
}

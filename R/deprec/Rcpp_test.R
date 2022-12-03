library("Rcpp")
cppFunction("bool isOddCpp(int num = 10) {
   bool result = (num % 2 == 1);
   return result;
}")
isoddCpp(42L)

# ax^2 + bx + c = 0

quad <- function(a,b,c){
  arg <- b^2 - 4*a*c
  if(arg >= 0){
  x1 <- (-b - sqrt(arg))/(2*a)
  x2 <- (-b + sqrt(arg))/(2*a)
  }
  else{stop("b^2 - 4ac < 0")}
  list(x1 = x1, x2 = x2)
}

quad(1,1,16)

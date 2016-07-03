library(lpSolve)

gurobi2lp <- function(model){
  res <- with(model, lpSolve::lp(direction=modelsense,
                                 objective.in = obj,
                                 const.mat = as.matrix(A),
                                 const.dir = sense,
                                 const.rhs = rhs
  ))
  result <- list(
    status = ifelse(res$status==0, 'OPTIMAL', 'other'),
    objval = res$objval,
    x = res$solution
  )
  return(result)
}
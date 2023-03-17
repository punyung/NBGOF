function (x, dist = "dnorm", fixed = NULL, link = NULL, start = NULL, 
          lower = NULL, upper = NULL, optimizer = "nlminb", control = NULL, 
          StdE_method = c("optim", "numDeriv"), silent = FALSE, ...) 
{
  #silent参数,负责warning的输出
  if (silent) 
    options(warn = -1)
  call <- match.call()
  #arguments里记录了分布函数，比如QNBD，的四个参数的名字
  arguments <- formals(dist)
  
  #control参数，暂时不知道有啥用，但如果不适list会报错
  if (!is.list(control)) {
    if (!is.null(control)) {
      stop("control argument must be a list \n \n")
    }
  }
  
  
  #dist参数，同上（这玩意是概率密度函数的名字）
  if (is.null(dist)) 
    stop("Distribution not specified \n \n")
  if (!is.character(dist)) 
    stop(paste0("'dist' argument must be a character ", 
                "string \n \n"))
  
  #产生了一个str向量solvers
  solvers <- unique(c("nlminb", "optim", "DEoptim", "ga", 
                      optimizer))
  
  
  #把optimizer赋值给了solvers
  solvers <- match.arg(optimizer, solvers)
  
  #处理link函数的输入参数，可以没有
  if (!is.null(link)) {
    if (length(match(link$over, names(arguments))) == 0) 
      #特定的分布要对应特定的link函数
      stop(paste0("Name(s) of linked parameter(s) do not agree with ", 
                  "arguments of ", dist, ". \n Please, change name(s) ", 
                  "specified in the entry 'over' of 'link' argument in \n", 
                  " function maxlogL.\n"))
    #over里存着parameters to map
    #fun里存着link function
    if (is.null(link$over) & !is.null(link$fun)) {
      warn <- paste0("You do not specify parameters to map, ", 
                     "however, you specify a link function \n ", 
                     "(the entry 'over' in 'link' argument is NULL ", 
                     "but the entry 'fun' is not NULL).\n")
      warning(warn)
    }
    
    if (!is.null(link$over) & is.null(link$fun)) 
      stop(paste0("You specify parameters to map, ", "however, you do not specify a link function \n", 
                  "(the entry 'fun' in 'link' argument is NULL ", 
                  "but the entry 'over' is not NULL).\n "))
  }
  #最好都给

  
  #和fixed相关的
    if (!is.null(fixed)) {
    if (length(match(names(fixed), names(arguments))) == 
        0) 
      stop(paste0("Name(s) of fixed (known) parameter(s) do not agree with ", 
                  "arguments of ", dist, ". \n Please, change names ", 
                  "specified in argument 'fixed' in function ", 
                  "maxlogL", "\n"))
    }
  
  
  #和输入的矩阵x相关的
  if (length(x) == 0 | is.null(x)) {
    stop(paste0("Vector of data is needed to perform maximum likelihood ", 
                "estimation. \n Please, specify the vector x in maxlogL ", 
                "function. \n"))
  }
  
  
  #把QNBD四个参数的名字输出为向量
  names_arguments <- names(arguments)
  
  pos_ncp <- sapply(names_arguments, function(y) grep("^ncp*.", 
                                                      y)[1])
  pos_ncp <- which(!is.na(pos_ncp))
  
  
  ###这里面的东西对于QNBD而言是不需要看的
  if (length(pos_ncp) > 0) {
    class_arguments <- sapply(arguments, class)
    num_ncp <- which((class_arguments[pos_ncp] == "numeric" | 
                        class_arguments[pos_ncp] == "symbol"))
    if (length(num_ncp) > 0) {
      fixed[[names_arguments[pos_ncp]]] <- arguments[[pos_ncp]]
    }
  }
  
  
  ###如果存在一个固定的参数
  names_fixed <- names(fixed)
  pos.deletion <- match(names_fixed, names_arguments)
  
  if (length(pos.deletion) > 0) 
    arguments <- arguments[-pos.deletion]
  nnum <- sapply(1:length(arguments), FUN = function(x) is.numeric(arguments[[x]])) # 判断参数是否为numeric
  nsym <- sapply(1:length(arguments), FUN = function(x) is.symbol(arguments[[x]])) # 判断参数是否为symbol
  npar <- length(nnum[nnum == TRUE]) + length(nsym[nsym == TRUE]) - 1 # 参数为numeric和symbol的数量相加再减1
  ll <- minus_lL(x = x, dist, dist_args = arguments, over = link$over, 
                 link = link$fun, npar = npar, fixed = fixed)
  
  
  
  if (is.null(lower))  # 设定下界
    lower <- rep(x = -Inf, times = npar)
  if (is.null(upper)) # 设定上界
    upper <- rep(x = Inf, times = npar)
  if (is.null(start))  # 开始运算的初始值
    start <- rep(x = 0, times = npar)
  if (!is.null(lower) & !is.null(upper) & !is.null(start)) {
    start <- link_apply(values = start, over = link$over, 
                        dist_args = arguments, npar = npar, link_fun = link$fun)
  }
  fit <- NULL
  if (optimizer == "nlminb") {
    nlminbcontrol <- control
    nlminb_fit <- nlminb(start = start, objective = ll, 
                         lower = lower, upper = upper, control = nlminbcontrol, 
                         ...)
    fit$par <- nlminb_fit$par
    fit$objective <- -nlminb_fit$objective
  }
  else if (optimizer == "optim") {
    optimcontrol <- control
    if (npar < 2) 
      optim_fit <- optim(par = start, fn = ll, lower = lower, 
                         upper = upper)
    optim_fit <- optim(par = start, fn = ll, control = optimcontrol, 
                       ...)
    fit$par <- optim_fit$par
    fit$objective <- -optim_fit$value
  }
  else if (optimizer == "DEoptim") {
    if (is.null(lower) | is.null(upper)) 
      stop("'lower' and 'upper'\n                                               limits must be defined\n                                               for 'DEoptim' optimizer", 
           "\n\n")
    DEoptimcontrol <- c(control, trace = FALSE)
    trace_arg <- which(names(DEoptimcontrol) == "trace")
    if (length(trace_arg) > 1) {
      control_index <- switch(as.character(call$control[[1]]), 
                              list = trace_arg[2])
      DEoptimcontrol[[control_index]] <- NULL
    }
    DE_fit <- DEoptim(fn = ll, lower = lower, upper = upper, 
                      control = DEoptimcontrol, ...)
    fit$original_fit <- DE_fit
    fit$par <- as.numeric(DE_fit$optim$bestmem)
    fit$objective <- -DE_fit$optim$bestval
  }
  else if (optimizer == "ga") {
    if (is.null(lower) | is.null(upper)) 
      stop("'lower' and 'upper'\n                                               limits must be defined\n                                               for 'GA::ga' optimizer", 
           "\n\n")
    plusll <- function(param) -ll(param)
    dots <- substitute(...())
    dots <- c(monitor = FALSE, control, dots)
    trace_arg <- which(names(dots) == "monitor")
    if (length(trace_arg) > 1) {
      dots$monitor <- NULL
    }
    ga_fit <- do.call(optimizer, c(list(type = "real-valued", 
                                        fitness = plusll, lower = lower, upper = upper), 
                                   dots))
    fit$original_fit <- ga_fit
    fit$par <- as.numeric(ga_fit@solution)
    fit$objective <- ga_fit@fitnessValue
  }
  else {
    fit <- do.call(optimizer, c(list(fn, lower, upper, start, 
                                     control, ...)))
  }
  fit$par <- link_apply(values = fit$par, over = link$over, 
                        dist_args = arguments, npar = npar, link_fun = link$fun)
  StdE_method <- match.arg(StdE_method, c("optim", "numDeriv"))
  ll.noLink <- try(minus_lL(x = x, dist, dist_args = arguments, 
                            over = NULL, link = NULL, npar = npar, fixed = fixed), 
                   silent = TRUE)
  if (StdE_method == "optim") {
    fit$hessian <- try(optim(par = fit$par, fn = ll.noLink, 
                             method = "L-BFGS-B", lower = fit$par - 0.5 * fit$par, 
                             upper = fit$par + 0.5 * fit$par, hessian = TRUE)$hessian, 
                       silent = TRUE)
    StdE_computation <- "Hessian from optim"
  }
  if (any(is.na(fit$hessian) | is.error(fit$hessian)) | any(is.character(fit$hessian)) | 
      StdE_method == "numDeriv") {
    fit$hessian <- try(numDeriv::hessian(ll.noLink, fit$par), 
                       silent = TRUE)
    StdE_computation <- "numDeriv::hessian"
  }
  if (any(is.na(fit$hessian) | is.error(fit$hessian)) | any(is.character(fit$hessian))) {
    StdE_computation <- paste0("'", StdE_method, "' failed")
    fit$hessian <- NA
    fit$StdE <- NA
  }
  else {
    fit$StdE <- sqrt(diag(solve(fit$hessian)))
  }
  names_numeric <- rep("", times = npar)
  j <- 1
  for (i in 1:length(arguments)) {
    if (is.numeric(arguments[[i]]) || is.symbol(arguments[[i]])) {
      names_numeric[j] <- names(arguments[i])
      j <- j + 1
    }
  }
  names_numeric <- names_numeric[-which(names_numeric == "x")]
  names(fit$par) <- names_numeric
  inputs <- list(call = call, dist = dist, fixed = fixed, 
                 link = link, optimizer = optimizer, start = start, lower = lower, 
                 upper = upper, data = x)
  outputs <- list(npar = npar - length(fixed), n = length(x), 
                  StdE_Method = StdE_computation, type = "maxlogL", par_names = names_numeric)
  result <- list(fit = fit, inputs = inputs, outputs = outputs)
  class(result) <- "maxlogL"
  if (silent) 
    options(warn = 0)
  return(result)
}

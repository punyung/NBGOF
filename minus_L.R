function (x, dist, dist_args, over, link, npar, fixed) 
{
# 输出：负对数似然函数值
# x: 表达矩阵
# dist: 概率密度函数
# dist_args: 其他分布函数参数 
# over: 需要进行约束的参数名称
# link：链接函数
# npar: 分布函数参数数量
# fixed: 固定的分布函数参数
  
  
  f <- function(param) {
    if (!is.null(link) & !is.null(over)) {
      linked_params <- link_list(over = over, dist_args = dist_args, 
                                 npar = npar)
      
      if (!is.null(linked_params)) {
        link_eval <- vector(mode = "list", length = length(linked_params))
        link <- paste0(link, "()")
        link_eval <- lapply(1:length(linked_params), 
                            FUN = function(x) eval(parse(text = link[x])))
        for (i in 1:length(linked_params)) {
          g_inv <- paste0("link_eval[[", i, "]]$g_inv")
          g_inv <- eval(parse(text = g_inv))
          param[[linked_params[i]]] <- do.call(what = "g_inv", 
                                               args = list(eta = param[[linked_params[i]]]))
        }
      }
    }
    logf <- do.call(what = dist, args = c(list(x = x), param, 
                                          log = TRUE, fixed))
    return(-sum(logf))
  }
  return(f)
}

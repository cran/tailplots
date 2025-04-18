###################################################################
## Paper: A Pareto Tail Plot Without Moment Restrictions by Klar ##
###################################################################



#This function returns t.hat evaluated at the sample points (and the ordered sample) mentioned on page 5 in the paper.
t_rec = function(x) {
    # recursive version, evaluated at sample points, starting with k=n-1
    x.vec = sort(x)
    n = length(x.vec)
    t = rep(0, n-1)
    s = rep(0, n-1)
    
    t[n-1] = (x.vec[n]-x.vec[n-1]) / (x.vec[n]+x.vec[n-1])
    
    x.j = x.vec[c(n-1,n)]
    for (k in (n-1):2) {
        # next smallest value
        x.m = x.vec[k-1] 
        t[k-1] = (n-k)/(n-k+2) * ( t[k] + 2/(n-k+1)/(n-k) * sum( (x.j-x.m) / (x.m+x.j)) )
        # add x.m to sample
        x.j = c(x.m,x.j) 
    }
    return( cbind(x.vec[-n], t) )
}

# This function is an efficient version of the function above with an (internal) double sum.
t_scalar = function(x, u) { 
    x.u = sort( x[x>=u] )
    n.u = length(x.u)
    s = 0
    x.j = x.u
    for (i in 1:(n.u-1))  {
        x.j = x.j[-1]
        s = s + sum( (x.j-x.u[i])/(x.u[i]+x.j) )
    }
    return( 2*s / (n.u*(n.u-1)) )
}

# This function applies a vector to t.scalar. 
t_vec = function(x, u) sapply(u, t_scalar, x=x)

# This function returns t.hat and both U-statistics mentioned on page 5 of the paper. 
t_comp = function(x, u) {
    n = length(x)
    x.u = sort( x[x>=u] )
    n.u = length(x.u)
    s = 0
    x.j = x.u
    for (i in 1:(n.u-1))  {
        x.j = x.j[-1]
        s = s + sum( (x.j-x.u[i])/(x.u[i]+x.j) )
    }
    return( c( s/choose(n.u,2), s/choose(n,2), choose(n.u,2)/choose(n,2) ) )
}


# This function computes 3 estimates for variance of t.hat. See section 4.1 of the paper for more details. 
var_pareto = function(x, u, method, R) { 
    x = sort(x)
    n = length(x)
    tco = t_comp(x,u)
    t = tco[1]
    u1 = tco[2]
    u2 = tco[3]

    if (method == "unbiased") {
        c1aa = c1ab = c1bb = 0
        c2aa = c2ab = c2bb = 0
        x.i = x
        #compute C_1^2 and C_2^2 in formula (9)
        for (i in 1:n)  { 
        x.i = x[-i]
        c1aa = c1aa + ( sum( (pmin(x[i],x.i) > u) * abs(x.i-x[i])/(x[i]+x.i) ) )^2
        c1ab = c1ab + sum( (pmin(x[i],x.i) > u) * abs(x.i-x[i])/(x[i]+x.i) ) * sum( (pmin(x[i],x.i) > u) )
        c1bb = c1bb + ( sum( (pmin(x[i],x.i) > u) ) )^2  
        
        c2aa = c2aa + sum( (pmin(x[i],x.i) > u) * ( abs(x.i-x[i])/(x[i]+x.i) )^2 )
        c2ab = c2ab + sum( (pmin(x[i],x.i) > u) * abs(x.i-x[i])/(x[i]+x.i) )
        c2bb = c2bb + sum( (pmin(x[i],x.i) > u) )  
        }
        #unbiased estimator
        var1 = (4*c1aa - 2*c2aa) / (n*(n-1)*(n-2)*(n-3)) - (4*n-6) * u1^2 / ((n-2)*(n-3)) 
        cov = (4*c1ab - 2*c2ab) / (n*(n-1)*(n-2)*(n-3)) - (4*n-6) * u1*u2 / ((n-2)*(n-3))
        var2 = (4*c1bb - 2*c2bb) / (n*(n-1)*(n-2)*(n-3)) - (4*n-6) * u2^2 / ((n-2)*(n-3)) 
        v.u = n / u2 * (var1 - 2*t*cov + t^2*var2)
        v.u = max(v.u, 0)   
        return( v.u )
    }

    if (method == "bootstrap") {
        boot = resample::bootstrap(x, t_scalar, R, args.stat = list(u))
        values = boot[[2]]
        #extract numeric values
        values = values[!is.na(values)] 
        v.b = n * u2 * var( values )
        return( v.b )
    }
    
    if (method == "jackknife") {
        jack = resample::jackknife(x, t_scalar, args.stat = list(u))
        values = jack[[2]]
        #extract numeric values
        values = values[!is.na(values)]
        v.j = (n-1)^2 * u2 * var( values )
        return( v.j )
    }
}

#This function is an optimitzed version of var.pareto. 
var_pareto_optimized = function(x, u, method, R) { 
    x = sort(x)
    n = length(x)
    tco = t_comp(x, u)
    t = tco[1]
    u1 = tco[2]
    u2 = tco[3]
    
    if (method == "unbiased") {
        # Precompute terms
        x_mat = matrix(rep(x, each = n), nrow = n)
        # Pairwise differences
        x_diff = abs(x_mat - t(x_mat))
        # Pairwise sums 
        x_sum = x_mat + t(x_mat)      
        pmin_mat = pmin(x_mat, t(x_mat))
        
        # Indicator matrix for condition pmin(x[i], x[j]) > u
        indicator = (pmin_mat > u) & lower.tri(pmin_mat, diag = FALSE)
        
        # Compute terms in a vectorized manner
        c1aa = sum(rowSums(indicator * x_diff / x_sum)^2)
        c1ab = sum(rowSums(indicator * x_diff / x_sum) * rowSums(indicator))
        c1bb = sum(rowSums(indicator)^2)
        
        c2aa = sum(rowSums(indicator * (x_diff / x_sum)^2))
        c2ab = sum(rowSums(indicator * x_diff / x_sum))
        c2bb = sum(rowSums(indicator))
        
        # Compute unbiased estimators
        denom = n * (n - 1) * (n - 2) * (n - 3)
        var1 = (4 * c1aa - 2 * c2aa) / denom - (4 * n - 6) * u1^2 / ((n - 2) * (n - 3))
        cov = (4 * c1ab - 2 * c2ab) / denom - (4 * n - 6) * u1 * u2 / ((n - 2) * (n - 3))
        var2 = (4 * c1bb - 2 * c2bb) / denom - (4 * n - 6) * u2^2 / ((n - 2) * (n - 3))
        
        # Final variance computation
        v.u = n / u2 * (var1 - 2 * t * cov + t^2 * var2)
        v.u = max(v.u, 0)
        return(v.u)
  }
  
    if (method == "bootstrap") {
        boot = resample::bootstrap(x, t_scalar, R, args.stat = list(u))
        values = boot[[2]]
        #extract numeric values
        values = values[!is.na(values)] 
        v.b = n * u2 * var( values )
        return( v.b )
  }
  
    if (method == "jackknife") {
        jack = resample::jackknife(x, t_scalar, args.stat = list(u))
        values = jack[[2]]
        #extract numeric values
        values = values[!is.na(values)] 
        v.j = (n-1)^2 * u2 * var( values )
        return( v.j )
  }
}

# This function computes the confidence intervals. See section 4.1 for more details. 
ci = function(x, u.vec, method, R=1000, conf.level = 0.95) {
    n = length(x)
    alpha = 1 - conf.level
    l.b = u.b = rep(0, length(u.vec)) 
    i = 1
    for (u in u.vec) {
        tco = t_comp(x,u)
        t.hat = tco[1]
        u2 = tco[3]
        vs = var_pareto_optimized(x, u, method, R)
        z = qnorm(1 - alpha/2)  
        l.b[i] = max( t.hat - z * sqrt(vs) / sqrt(n*u2), 0)
        u.b[i] = min( t.hat + z * sqrt(vs) / sqrt(n*u2), 1)
        i = i+1
    }
    return( rbind(l.b, u.b) )
}
#################
# Helpfunctions #
#################

# This function returs the theoretical value of t for the pareto distribution with shape a.  
t_pareto = function(a) {
    ifelse(a<=0, 1, a * (digamma((a+1)/2) - digamma(a/2)) - 1) 
  
}

#This function is used for root search. 
alpha_root = function(t, alpha.max=100) {
    f = function(u, t) t_pareto(u) - t 
    uniroot(f, c(0, alpha.max), t=t)$root
}
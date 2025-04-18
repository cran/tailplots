############################################################################
## Paper: A gamma tail statistic and its asymptotics by Iwashita and Klar ##
############################################################################



# This function returns g.hat and both U-statistics mentioned in equation (3) of the paper. 
g_comp = function(x, d) {
    x = sort(x)
    n = length(x)
    u1 = u2 = 0 
    i = 1; j = n
    while (i < j) {
        if (x[i]+x[j] > d) {
            u1 = u1 + sum( (x[j] - x[i:(j-1)]) / (x[j] + x[i:(j-1)])  )
            u2 = u2 + j - i
            j = j - 1
        } else {i = i + 1}
    }
    return( c(u1/u2, u1/choose(n,2), u2/choose(n,2)) )
}

# This function returns g.hat mentioned in equation (3) of the paper. 
g_scalar = function(x, d) {
    x = sort(x)
    n = length(x)
    u1 = u2 = 0 
    i = 1; j = n
    while (i < j) {
        if (x[i]+x[j] > d) {
            u1 = u1 + sum( (x[j] - x[i:(j-1)]) / (x[j] + x[i:(j-1)])  )
            u2 = u2 + j - i
            j = j - 1
        } else {i = i + 1}
    }
    return( u1 / u2 )
}

# This function applies a vector to g.scalar. 
g_vec = function(x,d) sapply(d, g_scalar, x=x)


# This function computes 3 estimates for variance of g.hat. See section 4 of the paper for more details. 
var_gamma = function(x, d, method, R) { 
    x = sort(x)
    n = length(x)
    gco = g_comp(x,d)
    g = gco[1]
    u1 = gco[2]
    u2 = gco[3]

    if (method == "unbiased") {
        c1aa = c1ab = c1bb = 0
        c2aa = c2ab = c2bb = 0
        x.i = x

        #compute C_1^2 and C_2^2 in formula (9)
        for (i in 1:n)  { 
            x.i = x[-i]
            c1aa = c1aa + ( sum( (x[i]+x.i > d) * abs(x.i-x[i])/(x[i]+x.i) ) )^2
            c1ab = c1ab + sum( (x[i]+x.i > d) * abs(x.i-x[i])/(x[i]+x.i) ) * sum( (x[i]+x.i > d) )
            c1bb = c1bb + ( sum( (x[i]+x.i > d) ) )^2  
            
            c2aa = c2aa + sum( (x[i]+x.i > d) * ( abs(x.i-x[i])/(x[i]+x.i) )^2 )
            c2ab = c2ab + sum( (x[i]+x.i > d) * abs(x.i-x[i])/(x[i]+x.i) )
            c2bb = c2bb + sum( (x[i]+x.i > d) )  
        }
        var1 = (4*c1aa - 2*c2aa) / (n*(n-1)*(n-2)*(n-3)) - (4*n-6) * u1^2 / ((n-2)*(n-3)) #unbiased estimator
        cov = (4*c1ab - 2*c2ab) / (n*(n-1)*(n-2)*(n-3)) - (4*n-6) * u1*u2 / ((n-2)*(n-3))
        var2 = (4*c1bb - 2*c2bb) / (n*(n-1)*(n-2)*(n-3)) - (4*n-6) * u2^2 / ((n-2)*(n-3)) 
        v.u = n / u2 * (var1 - 2*g*cov + g^2*var2)
        return( v.u )
    }
    
    if (method == "bootstrap") {
        boot = resample::bootstrap(x, g_scalar, R, args.stat = list(d))
        values = boot[[2]]

        #extract numeric values
        values = values[!is.na(values)] 
        v.b = n * u2 * var( values )
        return( v.b )
    }

    if (method == "jackknife") {
        jack = resample::jackknife(x, g_scalar, args.stat = list(d))
        values = jack[[2]]

        #extract numeric values
        values = values[!is.na(values)] 
        v.j = (n-1)^2 * u2 * var( values )
        return( v.j )
    }
}

# This function computes the confidence intervals. See section 5.2 for more details. 
ci_gamma = function(x, d.vec, method, R = 1000, conf.level = 0.95) {
    n = length(x)
    alpha = 1 - conf.level
    l.b = u.b = rep(0, length(d.vec)) 
    i = 1
    for (d in d.vec) {
        gco = g_comp(x,d)
        g.hat = gco[1]
        u2 = gco[3]
        vs = var_gamma(x, d, method, R)
        z = qnorm(1 - alpha/2)  
        l.b[i] = max( g.hat - z * sqrt(vs) / sqrt(n*u2), 0)
        u.b[i] = min( g.hat + z * sqrt(vs) / sqrt(n*u2), 1)
        i = i+1
    }
    return( rbind(l.b, u.b) )
}


#################
# Helpfunctions #
#################

# This function returs the theoretical value of g for the gamma distribution with shape a. See Remark 1 of the paper. 
g_gamma = function(a) {
    ifelse(a<=0, 1, ( 2^{2*a-1}*a*beta(a,a) )^{-1} )
}

#This function is used for the root search. 
alpha_root_gamma = function(g, alpha.max = 100) {
    f = function(a, g) g_gamma(a) - g
    uniroot(f, c(0, alpha.max), g=g)$root
}

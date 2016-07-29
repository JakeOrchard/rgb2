/*This mata file generates random data according to the specified distribution.

Author --Jacob Orchard

v1.2
Update 7/29/2016*/



version 13
mata:
	function rgb2(string myvar, scalar delta, scalar sigma, scalar p, scalar q)
	{
	nobs = st_nobs()
	base = runiform(nobs,1)
	newvar = J(nobs,1,0)
	paravec = delta,sigma,p,q
	
	if (nobs <1000){
	
		for (i=1; i <=nobs; i++){
			newvar[i] = estimategb2(paravec,base[i])
			}
		}
	else{
		cdflist = interpolategb2(paravec)
		
		for (i=1; i <=nobs; i++){
			newvar[i] = nearest_gb2(cdflist,base[i])
			}
	
	}
	
	(void) st_addvar("double", myvar)
	st_store(.,myvar,newvar[.,1])
	}
	end
	
mata:
	function estimategb2(paravec,unifval)
	{
	
		delta = paravec[1]
		sigma = paravec[2]
		p     = paravec[3]
		q     = paravec[4]
		
		S = optimize_init()
		optimize_init_which(S,"min")
		optimize_init_technique(S,"dfp nr")
		optimize_init_evaluator(S,&closestgb2())
		
		//Expected Value (starting value)
		start = exp(delta)*( (exp(lngamma(p+sigma))*exp(lngamma(q-sigma)))/( exp(lngamma(p))*exp(lngamma(q))))
		optimize_init_params(S,start)
		optimize_init_argument(S,1,paravec)
		optimize_init_argument(S,2,unifval)
		_optimize(S)
	    ans    = optimize_result_params(S)
		return(ans)
	
	}
	
end

mata:
	function interpolategb2(paravec)
	{
		cdflist = J(1000,1,0)
		
		for (i=1; i <=1000; i++){
		cdflist[i] = estimategb2(paravec,(i/1000))
		}
	return(cdflist)
    }
end

mata:
	function nearest_gb2(cdflist,unifval)
	{
		nx = 1000*unifval
		intx = floor(nx)
		remx = nx - intx
		
		closegb2 = cdflist[intx]*(1-remx) + cdflist[intx+1]*remx
		return(closegb2)
	}
end

mata: 
	function closestgb2(real scalar todo, real rowvector p, paravec,unifval, v,g,H)
	{
	
	v = abs(gb2_cdf(p,paravec)-unifval)
	
	}
	end
	
mata:
function gb2_cdf(matrix y, paravec)
	{
	ones = J(1,cols(y),1)
	delta = paravec[1]
	sigma = paravec[2]
	p = paravec[3]
	q = paravec[4]
	zu = ((y:/exp(delta)):^(1/sigma)):/(ones+(y:/exp(delta)):^(1/sigma))
	F =  ibeta(p,q,zu)
	return(F)
	}
end
	
	

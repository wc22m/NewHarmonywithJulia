# define the harmony function
using Statistics
function harmony(x)
	x/(1.1+x)
	end
function harmonyInverse(h)
# h= x/(1.1+x)
# 1.1h +hx =x
# 1.1h = x- hx
# 1.1h = (1-h)x
# 1.1h/(1-h) =x
  return (1.1*h)/(1-h)
end
mincoeff=0.01

maxiter=500
verbose = false
si=[100.0,90,80,70,60]# initial capital stock
labourproductivity=[2.173913043, 2.19,2.21,2.23,2.25]
labouravailable=[25,25.125,25.250625,25.37687813,25.50376252]
g=[33,33.165,33.330825,33.49747913,33.66496652]
labourproductivity=[2.173913043
	2.173913043
	2.273913043
	2.313913043
	2.343913043
]
investments=[0.0,0.0,0.0,0.0,0.0]
capitalproductivity=
[
	0.5
	0.5
	0.5
	0.5
	0.5
]
function computeOutput(capital,labour)
	return min(capital .* capitalproductivity,labour .* labourproductivity)
end
o= computeOutput(si,labouravailable)
fulfill=o ./ g
h= [harmony(x) for x in fulfill]
println("initial harmony without capital replacement",h)
depreciationhorizon = 14
epsilon = 3/depreciationhorizon
# we use a linear depreciation to find out how much capital we need some years earlier 
function inversedepreciate(capital, yearsearlier)
   if yearsearlier ==1 
      return capital
   end
   return capital * depreciationhorizon/(depreciationhorizon -(yearsearlier -1))
end
function depreciateamountbyyears(amount,years)
	if years >depreciationhorizon
		return 0
	end
	return amount * (depreciationhorizon - years )/depreciationhorizon
end
function computeMinGrossOutputForGoal(goal)
	ming = goal
	for i= 1:15
		mincapreq= ming ./ capitalproductivity
		depreciationallowance = mincapreq ./ depreciationhorizon
		ming = goal + depreciationallowance
	end
	return ming
end
TheLastYear=5
 #=
 
1 For the last year of the plan set a net output target such that \label{item:settarget}
 
    a Gross output is such as to ensure full employment of the workforce
    b Sufficient investment is being carried out to compensate for depreciation during the year.
  
 This net output target must be a scaled version of the originally specified
 target. If the original target would have caused unemployment it will be up-scaled
 otherwise down-scaled.
 How do we compute the full employment net output  which we will designate as $N$ here?
 We are assumed to know the productivity $k$ of ( each kind of ) fixed capital,
 the depreciation rate of fixed capital $\Delta$ and the productivity of labour $a$.
 We have the relationships for the gross output $O$ as follows
 $$O\le kC$$ $$ O \le aL $$ $$N=O-\Delta C$$ where $C$ is the fixed capital stock and $L$ the labour used.
 So \begin{equation}N\le kC -\Delta C
      \label{eq:nle}
 \end{equation} $$\Delta C \le kC -N$$ $$N\le aL -kC -N$$
 so $$kC \le aL$$ let's assume equality and allow $C = \frac{a}{k}L$ then by
 (\ref{eq:nle}) we have $$N=(k-\Delta)\frac{a}{k}L$$
 =#
 function For_the_last_year_of_the_plan_set_a_net_output_target(last_year)
   k=capitalproductivity[last_year]
   Delta = 1.0 /depreciationhorizon
   a=labourproductivity[last_year]
   L=labouravailable[last_year]
   N=(k-Delta)*(a/k)*L
   return N
 end
 g[TheLastYear]=For_the_last_year_of_the_plan_set_a_net_output_target(TheLastYear)
 #=
2 Assign to each year's fixed capital stock the starting stock depreciated by 
 the appropriate depreciation rate on that type of fixed capital stock.
 =#
 function Assign_to_each_year_capital_stock()
	startingstock=si[1]
	stock = 
    for year = 1:TheLastYear
     if (year -1 )<depreciationhorizon
       global si[year] = si[1] *(depreciationhorizon -(year -1))/depreciationhorizon
     else
		global si[year]=0
     end
    end
 end
 
 Assign_to_each_year_capital_stock()
 
 println("Depreciated initial capital stocks ",si)
 println("For the last year of the plan set a net output target ",For_the_last_year_of_the_plan_set_a_net_output_target(TheLastYear),g)
 #=
 3 Compute the degree of goal fulfillment that is in principle possible
 given these fixed capital stocks and the available workforce. From this compute
 the Harmonies.
 \label{item:goalfulfill}
 =#
 doagain =true
 count=1
 while doagain
	 function Compute_the_degree_of_goal_fulfillment_that_is_in_principle_possible(s,l,g )
	    computeOutput(s,l)./g
	 end
	 # note that in_principle_possible is degree of goal fulfillment not level of output
	 in_principle_possible=Compute_the_degree_of_goal_fulfillment_that_is_in_principle_possible(si,labouravailable,g )
	 
	 if verbose
		println("Goal fulfillment in principle posssible ",in_principle_possible)
	 end
	 
	 #=
	4 Compute the mean harmony, and standard deviation of the harmonies over the whole time period.
	 
	 =#
	 function computeHarmonies()
		 [harmony(x) for x in ( in_principle_possible-(investments ./ g))]
		
	 end
	 global h =computeHarmonies()
	 function Compute_mean_harmony_and_standard_deviation()
	     harmonies = computeHarmonies()
	     return(mean(harmonies),std(harmonies))
	 end
	  meanh,stdh = Compute_mean_harmony_and_standard_deviation()
	  
	  if verbose
		 println("Harmonies ",computeHarmonies()) 
		 println(" mean and standard deviation of harmonies ",Compute_mean_harmony_and_standard_deviation()," c of v ",stdh/meanh)
	  end
	  
	  
	 #=
	5 If the coefficient of harmony variation has fallen below some threshold, terminate the plan prepartion, otherwise continue to   step 6.
	=#
	if (stdh/meanh < mincoeff)||(count>maxiter)
	        println("terminate after ",count," iterations with coefficient of variation = ",stdh/meanh)
	        
			println("final capital stocks ",si)
			println("final harmonies ",h)
			println("Goal fulfillment in principle posssible ",in_principle_possible)
			gross=computeOutput(si,labouravailable)
			println("gross output ",gross)
			employment=gross ./ labourproductivity
			println("labour usage ",employment )
			println("percentage employment ",100* employment ./ labouravailable)
			global doagain=false
	end
	global count =count+1
	if doagain 
		#=
		6 Select the year with the lowest harmony. 
		 =#
		 function Select_the_year_with_the_lowest_harmony()
		     H=computeHarmonies();
		     minim=1000
		     minyear=TheLastYear
		     for year = 1:TheLastYear
		       if minim > H[year]
		          minim=H[year]
		          minyear=year
		       end
		     end
		     return minyear
		 end
		 
		lowyear=Select_the_year_with_the_lowest_harmony()
		
		if verbose
		 println(" lowest harmony year ",lowyear)
		end
		 #=
		7 Estimate by how much the  production  of this year would have to be scaled up in order to be at the mean harmony.
		This is returned as a ratio.
		To do this we first work out the level of net production needed to meet the mean harmony. This involves inverting the harmony function.
		 =#
		 function Estimate_how_much__production_to_be_scaled_up(this_year)
		    targetgoalfulfillmentlevel = harmonyInverse(meanh)+investments[this_year]/g[this_year]
		    return targetgoalfulfillmentlevel / in_principle_possible[this_year]
		 end
		 
		 
		 upscale=Estimate_how_much__production_to_be_scaled_up(lowyear )
		 
		 if verbose
			println(" how much this year needs to be scaled by ", upscale)
		 end
		 #=
		8 Attempt to scale up some fraction $\epsilon$ of the way towards the mean harmony by \label{item:realloc}
		 scheduling investment in  at least one previous year. This will obviously alter both the consumption in previous years
		  and the fixed capital stocks in one or more years.
		  =#
		  function update_subsequent_years_capital(firstyearavailable,amount)
		      finalyear = TheLastYear
		      if finalyear <firstyearavailable+depreciationhorizon
		         finalyear=firstyearavailable+depreciationhorizon
		      end
		      for y= firstyearavailable:TheLastYear
		        si[y]=si[y]+depreciateamountbyyears(amount,y-firstyearavailable,)
		      end
		 end
		 function Attempt_to_scale_up(year,byfractionofoutput)
		   # determine how much additional capital we need to achieve this
		   additionalcapital = si[year] * (byfractionofoutput -1)*epsilon
		   # now find which of the previous years will supply this
		   gaininharmonyifachieved = harmony(in_principle_possible[year]*byfractionofoutput-investments[year]/g[year])-harmony(in_principle_possible[year]-investments[year]/g[year])
		   yearsPrevious=[y for y =1:year-1]
		   bestyear=0; bestgain = -100000000.0; lossinbestyear=0
		   for y in yearsPrevious
		     function lossfromInvesting(inyear,amountmore)
		         initialinvestmentfraction=investments[inyear]/g[inyear]
		         additionalfraction=amountmore/g[inyear]
		         currentHarmony = harmony(in_principle_possible[inyear]-initialinvestmentfraction)
		         futureHarmony= harmony(in_principle_possible[inyear]-initialinvestmentfraction-additionalfraction)
		         return   futureHarmony -currentHarmony
		     end
			 loss =lossfromInvesting(y,inversedepreciate(additionalcapital,year-y))# negative
			 if (loss+gaininharmonyifachieved)>bestgain
				bestyear=y
				lossinbestyear=loss
				bestgain = loss+gaininharmonyifachieved
			 end
		   end
		   if verbose
			   println(" capital wanted  for fraction ",1/depreciationhorizon," of this is ",additionalcapital," to achieve gain ",gaininharmonyifachieved)
			   println(" take from year ",bestyear," which will loose ",lossinbestyear)
			   println(" previous capital stocks ",si)
		   end
		   investments[bestyear]=investments[bestyear]+inversedepreciate(additionalcapital,year-bestyear)
		   update_subsequent_years_capital(bestyear+1,inversedepreciate(additionalcapital,year-bestyear))
		   if verbose 
			println(" updated capital stocks ",si)
		   end
		 end
		 
		Attempt_to_scale_up(lowyear,1+(upscale-1)/depreciationhorizon)
		  #=
		9 Go back to step \ref{item:goalfulfill}.
		  
		 
		 
		 =#
	 end
end

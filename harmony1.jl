# define the harmony function
function harmony(x)
	x/(1.1+x)
	end

si=[100,90,80,70,60]# initial capital stock
labourproductivity=[2.173913043, 2.173913043,2.173913043,2.173913043,2.173913043]
labouravailable=[25,25.125,25.250625,25.37687813,25.50376252]
g=[33,33.165,33.330825,33.49747913,33.66496652]
labourproductivity=[2.173913043
	2.173913043
	2.273913043
	2.313913043
	2.343913043
]
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
depreciationhorizon = 10
function computeMinGrossOutputForGoal(goal)
	ming = goal
	for i= 1:15
		mincapreq= ming ./ capitalproductivity
		depreciationallowance = mincapreq ./ depreciationhorizon
		ming = goal + depreciationallowance
	end
	return ming
end
mingross=computeMinGrossOutputForGoal(g) 

println("stable capital stock gross output to meet original goal",mingross)
sn = mingross ./ capitalproductivity 

println("si ",si)
println("sn ",sn)
o2= computeOutput(sn,labouravailable)
fulfill2 = o2 ./g
 
labuse=mingross ./labourproductivity
unemp = (labouravailable - labuse)./ labouravailable
println(" unemployment rate assuming sn",unemp)
# now compute the maximum output that could
# be achieved with full employment each year and 
# scale up the goals
fullemploymentgoal = g .*(1 .+unemp)
println("Original goal ",g)
println("Full employment goal ",fullemploymentgoal)
fullempgrossOutput=computeMinGrossOutputForGoal(fullemploymentgoal)
println("full emp gross output ", fullempgrossOutput)
println("labour and capital constrained output attempting to meet original goal",o2)
println("labour constrained max output ",computeOutput(1000 .* sn,labouravailable))

 

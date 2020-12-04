# file that reads in the csv files and runs plan
using DelimitedFiles
using LinearAlgebra
using Statistics
#-----------------------------------------------------------------------
# Constants
#-----------------------------------------------------------------------
mincoeff = 0.034
maxiter = 3000
verbose = false
linear = false
inversedepreciateinvestments = false
depreciationhorizon = 14 # hopefully this is a long enough estimate
initialinvestmentlevel = 0.7
epsilon = 0.25 / depreciationhorizon
#-----------------------------------------------------------------------
# type declarations
#-----------------------------------------------------------------------
struct PlanProblem
    # data as read in
    headers::Any         # headers from the labtargs csv file
    flows::Any# data from the flow io table csv file
    caps::Any# data from the capital stock  csv file
    dep::Any# data from the depreciation csv file
    labtarg::Any  # data from the labour supply and plan target file
    #  key matrices
    A::Any  # technology matrix= generateAMatrix(flows)
    V::Any  # labour values = compute_Labour_Value_From_Augmented_A_Matrix(A)
    C::Any  # capital matrix = generateCMatrix(flows, caps)
    D::Any  # depreciation matrix = generateDMatrix(C, dep)


    AD::Any # = A + D
    IminusA::Any     # = I - A
    InvIA::Any       # = inv(IminusA)
    #  targets and labour supply
    labouravailable::Any # labour available each year= (extractlaboursupply(labtarg, headers))
    g::Any # explicit plan targets by year = extracttargets(labtarg, headers)
    TheLastYear::Any # number of years to run plan = size(g)
    caprows::Any
    capcols::Any         # size of capital table
end
mutable struct Scenario
    prob::Any# the current plan problem
    O::Array{Float64,2} # gross output to meet final demand including investments at this stage of algorithm [year,product]
    si::Array{Float64,3} # the capital stock tensor at this stage of algorithm [year,srcprod,useprod]
    investments::Array{Float64,3} # the investment tensor at this stage of algorithm [year,srcprod,useprod]
    goal_fullfilment_ratio_vector::Any      # indexed by years, the share of O that can actually be met
    h::Array{Float64,1} # harmonies by year
    meanh::Any# mean harmony
    stdh::Any# harmony std deviation
    investmentsByTypeAndYear::Any# investments totaled by type of good [year,type]
    netoutputs::Any                          # output - intermediate consumption - investment
    targets::Any# final net output targets including investment
end
#-----------------------------------------------------------------------
# global variables
#-----------------------------------------------------------------------
doagain = true
#-----------------------------------------------------------------------
# Pure functions section
#-----------------------------------------------------------------------

function readinspreadsheets(flow, cap, dep, labtarg)
    flowsd, flowsh = DelimitedFiles.readdlm(flow, ',', Float64, header = true)
    capsd, capsh = readdlm(cap, ',', Float64, header = true)
    depsd, depsh = readdlm(dep, ',', Float64, header = true)
    labtargsd, labtargsh = readdlm(labtarg, ',', Float64, header = true)
    return labtargsh, flowsd, capsd, depsd, labtargsd
    # returns a tuple (labtargcolheaders,flowmatrix, capitalmatrix, depreciation matrix, targetmatrix )
end
function extractlaboursupply(labtargs, headers)
    hrow, hcol = size(headers)
    lrow, lcol = size(labtargs)
    if headers[1, hcol] != "Labour"
        throw(error("Labour not last header of labtargs"))
    end
    l = zeros(Float64, lrow + depreciationhorizon, 1)
    for i = 1:lrow+depreciationhorizon
        if i < lrow
            l[i] = labtargs[i, hcol]

        else
            l[i] = labtargs[lrow, hcol]
        end
    end
    return l
end
function extracttargets(labtargs, headers)
    hrow, hcol = size(headers)
    lrow, lcol = size(labtargs)
    if headers[1, 1] != "Year"
        throw(error("Year not first header of labtargs"))
    end
    v = @view labtargs[:, 2:(hcol-1)]
    g = zeros(Float64, lrow + depreciationhorizon, hcol - 2)
    g[1:lrow, :] = v[:, :]
    for i = lrow+1:lrow+depreciationhorizon
        g[i, :] = v[lrow, :]
    end
    return hcat(g, zeros(Float64, lrow + depreciationhorizon, 1))
end
# the D specifies how much depreciation must be allowed for given steady
# state production of a given output vector

function generateDMatrix(C, dep)
    crow, ccol = size(C)
    drow, dcol = size(dep)
    D = zeros(Float64, crow, ccol)
    D[1:drow, 1:dcol] = C[1:drow, 1:dcol] .* dep
    D
end
function generateCMatrix(flowmatrix, cap)
    # assume last row is output
    # assume second to last row is labour
    # first get the size
    rows, cols = size(flowmatrix)
    N = rows - 1
    if N != (cols + 1)

        throw(error("rows not equal to cols +2 in flowmatrix"))
    end
    # create a  matrix view of the outputinvestments
    O = @view flowmatrix[rows, 1:cols]
    #  create a normalised version of the rest
    Normcap = zeros(Float64, N, N)
    for row = 1:N-1
        for col = 1:N-1
            Normcap[row, col] = cap[row, col] / O[col]
        end
    end
    # return an N by N D matrix
    return Normcap
end
function generateAMatrix(flowmatrix)
    # assume last row is output
    # assume second to last row is labour
    # first get the size
    rows, cols = size(flowmatrix)
    N = rows - 1
    if N != (cols + 1)
        throw(error("rows not equal to cols +2 in flowmatrix"))
    end
    # create a  matrix view of the output
    O = @view flowmatrix[rows, 1:cols]
    #  create a normalised version of the rest
    Normfloe = zeros(Float64, N, cols)
    for row = 1:N
        for col = 1:cols
            Normfloe[row, col] = flowmatrix[row, col] / O[col]
        end
    end
    # return an N by N A matrix
    return hcat(Normfloe, zeros(Float64, N, 1))
end
function inversedepreciate(capital, yearsearlier, depreciationmatrix)
    if yearsearlier <= 1
        return capital
    end
    #println("inversedepreciate($capital , $yearsearlier , $depreciationhorizon ")
    c = inversedepreciate(capital, yearsearlier - 1, depreciationmatrix)
    return c ./ (1.0 .- depreciationmatrix)
end

function rowiseScale(A, b)
    N, cols = size(A)
    C = zeros(Float64, N, cols)
    for row = 1:N
        for col = 1:cols
            C[row, col] = A[row, col] * b[col]
        end
    end
    return C
end
function grossOutputForDemandf(f::Array{Float64,1}, A::Array{Float64,2})
    invia = inv(I - A)
    gross = invia * f
    #display(invia)
    #println("\ndemand ");println(f);println("\n needs gross") ;println(gross )
    return gross
end;
function computeMaxPossOutput_in_principle_possible(s, l, A, C, O, year)
    if verbose
        print("\ncomputeMaxPossOutput  with stocks ")
        display(s[year, :, :])
    end
    N, n = size(A)
    # first get the labour constraint
    labrow = A[N, :]
    thisyearoutput = O[year, :]
    labourneeded = sum(labrow .* thisyearoutput)
    labourconstraint = l[year] / labourneeded
    if verbose
        println("\nC")
        display(C)
        println("\nthis year output", thisyearoutput)
    end
    # next the capital constraint
    capneeded = rowiseScale(C, thisyearoutput)
    if verbose
        println("\ncapneeded ")
        display(capneeded)
    end
    caprow, capcol = size(s[year, :, :])
    metratios = s[year, :, :] ./ capneeded[1:caprow, 1:capcol]


    capitalmin = min(metratios...)
    overallmin = min(capitalmin, labourconstraint)
    if verbose
        println("\ncapitalmetratio year ", year)
        display(metratios)
        println("labourconstraint ", labourconstraint)
        println(year, " ", overallmin)
    end

    return overallmin
end
function computeMaxPossOutputRatioAllYears(s, l, A, C, O)

    lastyear, cols = size(O)
    # println("lastyear cols", lastyear, cols)
    return [
        computeMaxPossOutput_in_principle_possible(s, l, A, C, O, year)
        for year = 1:lastyear
    ]
end
# note that in_principle_possible_goal_fullfilment_ratio_vector is degree of goal fulfillment not level of output
# s= capital stocks
# l labour available
# O final output target
# A technology matrix
# C capital matrix
function Compute_the_degree_of_goal_fulfillment_that_is_in_principle_possible_goal_fullfilment_ratio_vector(
    s,
    l,
    O,
    A,
    C,
)

    return computeMaxPossOutputRatioAllYears(s, l, A, C, O)
end
function compute_Labour_Value_From_Augmented_A_Matrix(A)
    # this returns a vector of labour values using an A matrix in which it is assumed
    # that the last row is labour input
    N, n = size(A)
    v = zeros(Float64, N)
    v[N] = 1
    for steps = 1:12
        nv = A' * v
        v = nv
        v[N] = 1
    end
    return v
end
function inproduct(v1, v2)
    v1 = vec(hcat(v1...))
    v2 = vec(hcat(v2...))
    vs1 = size(v1)
    vs2 = size(v2)
    short = min(vs1[1], vs2[1])
    t = 0.0
    for i = 1:short
        t = t + v1[i] * v2[i]
    end
    t
end
function valueof(valuation::Array{Float64,1}, commodityvector::Array{Float64,1})
    return inproduct(valuation, commodityvector)
end
function valueof(valuation::Array{Float64,1}, commoditymatrix::Array{Float64,2})
    n, m = size(commoditymatrix)
    subv = @view valuation[1:n]
    #println("commoditymatrix =");display(commoditymatrix);
    #println("sub valuation matrix =",subv)
    return sum(commoditymatrix' * subv)
end
function For_the_last_year_of_the_plan_return_a_net_output_target(last_year, A, D, g, l)
    N, n = size(A)
    # construct a modified A allowing for capital requirements


    Gross = grossOutputForDemandf(g[last_year, :], A)

    labourRow, cols = size(A)
    labourUsed = Gross' * A[labourRow, :]
    depreciationvec = D * Gross
    # scale gross output to bring about full employment
    #    println("labour used ", labourUsed)
    newtarget = (g[last_year, :] + depreciationvec) .* (l[last_year] / labourUsed)
    scaledgross = Gross * (l[last_year] / labourUsed)
    depreciationallowanceMatrix = rowiseScale(D, scaledgross)

    return (depreciationallowanceMatrix, newtarget)
end

function display_result(
    year::Int64,
    investmentsbytype::Array{Float64,2},
    C,
    si::Array{Float64,3},
    harmonies::Array{Float64,1},
    headings,
    gross::Array{Float64,1},
    A::Array{Float64,2},
    labsupply,
    goal,
    s,
)

    println("\n Input output table year,", year, ", harmony,", harmonies[year])
    for h in headings[2:end]
        print(",", h)
    end

    println(",IntermediateConsumption,Investment,FinalConsumption,ConsumptionTarget")
    rows, cols = size(A)
    labused = 0
    invyears, invrows = size(investmentsbytype)
    for row = 1:rows
        print(headings[row+1])
        tot = 0.0
        toti = 0
        for col = 1:cols
            cell = A[row, col] * gross[col]
            print(",", cell)
            tot = tot + cell
        end
        if row < invrows
            toti = investmentsbytype[year, row]
        end
        println(
            ",",
            tot,
            ",",
            toti,
            ",",
            (row < rows ? s.netoutputs[year][row] : ""),
            ",",
            (row < rows ? goal[year, row] : ""),
        )
        if row == rows
            labused = tot
        end
    end
    print("GrossOutput")
    for col = 1:cols
        print(",", gross[col])
    end
    println()
    println("Capital stock")
    rows, cols = size(C)

    for row = 1:rows
        print(headings[row+1])

        for col = 1:cols
            cell = C[row, col] * gross[col]
            print(",", cell)

        end
        println()
    end
    println("employment ratio,", labused / labsupply)
end
#=
  Assign to each year's fixed capital stock the starting stock depreciated by
 the appropriate depreciation rate on that type of fixed capital stock.
 =#

function Assign_to_each_year_capital_stock(caps, deps, TheLastYear)
    caprow, capcol = size(caps)

    si = zeros(Float64, TheLastYear, caprow, capcol)
    si[1, :, :] = caps[:, :]
    for i = 2:TheLastYear
        si[i, :, :] = depreciateamountbyyears(caps, i, deps)[:, :]
    end

    return si

end

function harmony(x)
    if linear
        return x
    end
    x / (1.1 + x)
end
function harmonyInverse(h)
    # h= x/(1.1+x)
    # 1.1h +hx =x
    # 1.1h = x- hx
    # 1.1h = (1-h)x
    # 1.1h/(1-h) =x
    if linear
        return h
    end
    y = (1.1 * h) / (1 - h)
    return y
end
function finaloutputforGross(gross::Array{Float64,1}, A::Array{Float64,2})
    return (I - A) * gross
end

function update_outputs(s::Scenario)
    lastyear = s.prob.TheLastYear
    investmentsbytype = s.investmentsByTypeAndYear = investmentsByTypeandYear(s.investments)

    s.targets = s.prob.g[:, 1:end-1] + s.investmentsByTypeAndYear
    gplusi = zeros(Float64, size(s.prob.g))
    gplusi[:, 1:end-1] = s.targets

    tmpO = [grossOutputForDemandf(gplusi[i, :], s.prob.A)' for i = 1:s.prob.TheLastYear]
    # make sure it is a matrix
    s.O = vcat(tmpO...)
    Compute_the_degree_of_goal_fulfillment_that_is_in_principle_possible_goal_fullfilment_ratio_vector(
        s,
    )
    possibleratios = s.goal_fullfilment_ratio_vector
    # println("possible ratios",possibleratios);
    #println("\n investments");display( investmentsbytype)
    #println("\n targets");display(s.targets)
    finaloutput = [s.targets[i, :] .* possibleratios[i] for i = 1:lastyear]

    # println("\nfinal outputs");display(finaloutput)
    netoutputs = [finaloutput[i] .- investmentsbytype[i] for i = 1:lastyear]
    s.netoutputs = netoutputs

    # println("\nnetoutputs");display(netoutputs)
end
function computeHarmonies(
    s::Scenario,
    possibleratios,
    investmentsbytype::Array{Float64,2},
    netgoals::Array{Float64,2},
    O,
    A,
)
    lastyear = s.prob.TheLastYear


    fulfillmentratio = [s.netoutputs[i] ./ s.prob.g[i, 1:end-1] for i = 1:lastyear]
    # println("\n fulfillment ratio")
    # display(fulfillmentratio)
    harmonyarrayofarrays = [map(x -> harmony(x), fulfillmentratio[i]) for i = 1:lastyear]
    #println("\n harmony matrix")
    #display(harmonyarrayofarrays)
    # the "end-1" below is because we ignore labour when evaluating harmony
    minimised = [min(((harmonyarrayofarrays[i])[1:end-1])...) for i = 1:lastyear]
    # println("\n minimised to ")
    # display(minimised)
    return minimised
end

function Compute_mean_harmony_and_standard_deviation(harmonies)
    return (mean(harmonies), std(harmonies))
end
function computeHarmonies(s::Scenario)
    s.h = computeHarmonies(
        s,
        s.goal_fullfilment_ratio_vector,
        s.investmentsByTypeAndYear,
        s.prob.g,
        s.O,
        s.prob.A,
    )
    meanh, stdh = Compute_mean_harmony_and_standard_deviation(s.h)
    s.meanh = meanh
    s.stdh = stdh
end

function depreciateamountbyyears(amount, years, depreciationmatrix)
    # depreiciate using exponential depreciation schedule
    if years <= 1
        return amount
    end
    return depreciateamountbyyears(amount, years - 1, depreciationmatrix) .*
           (1.0 .- depreciationmatrix)
end

# reduce a rank 3 tensor of investments to a matrix[year,sourcindustry]
function investmentsByTypeandYear(investments)
    years, sources, dests = size(investments)
    bytypeandyear = zeros(Float64, years, sources)
    for y = 1:years
        for s = 1:sources
            for d = 1:dests
                bytypeandyear[y, s] = bytypeandyear[y, s] + investments[y, s, d]
            end
        end
    end
    return bytypeandyear
end
function Select_the_year_with_the_lowest_harmony(H::Array{Float64,1})
    # we do not select any year for which the goal achievement <=0
    minim = 1000
    minyear = 1
    years, = size(H)
    for year = 1:years
        println("H[$year]=", H[year], " harmony inverse ", harmonyInverse(H[year]))
        if minim > H[year]# && harmonyInverse(H[year]) > 0
            minim = H[year]
            minyear = year
        end
    end
    #println("minyear ",minyear)
    return minyear
end
function update_subsequent_years_capital(firstyearavailable, amount, horizons, si)
    years, rows, cols = size(si)

    for y = firstyearavailable:years
        si[y, :, :] =
            si[y, :, :] + depreciateamountbyyears(amount, y - firstyearavailable, horizons)
    end
end
function setup_preliminary_investment_schedule(
    year,
    investments,
    preassignedcapital,
    horizons,
    si,
)
    investments[year, :, :] = preassignedcapital[:, :]
    update_subsequent_years_capital(year + 1, preassignedcapital, horizons, si)
end;

function updateScenario(s::Scenario, bestyear, capital)
    s.investments[bestyear, :, :] = s.investments[bestyear, :, :] + capital
    if sum(capital) == 0
        throw(error("attempt to update with zero investment "))
    end
    update_subsequent_years_capital(bestyear + 1, capital, s.prob.dep, s.si)
    if verbose
        println("\n bestyear $bestyear ")
        display(s.investments[bestyear, :, :])
        println(" updated capital stocks ", s.si)
    end
    return s
end
function variant_scenario(s::Scenario, newinvestmentforyear, year)
    s2 = deepcopy(s)
    updateScenario(s2, year, newinvestmentforyear)
    update_outputs(s2)
    computeHarmonies(s2)
    return s2
end
function gainfromInvesting(s::Scenario, fromyear::Integer, amountmore)
    variant = variant_scenario(s, amountmore, fromyear)
    futureHarmony = variant.meanh
    currentHarmony = s.meanh

    return (futureHarmony - currentHarmony, variant)
end
function capitalbyyear(s::Scenario)
    [sum(s.si[i, :, :]) for i = 1:s.prob.TheLastYear]
end
function Attempt_to_scale_up(
    s::Scenario,
    destyear::Integer,
    scaleincrementtomeetmeanharmony::Float64,
)
    # determine how much additional capital we need to achieve mean harmony
    V = s.prob.V
    g = s.prob.g
    csy = s.si[destyear, :, :]
    n, m = size(csy)
    additionalcapital = csy .* scaleincrementtomeetmeanharmony
    # remove negative values since some capital stocks may already be higher
    for i = 1:n
        for j = 1:m
            if additionalcapital[i, j] < 0
                additionalcapital[i, j] = 0
            end
        end
    end

    if verbose
        println("additional capital for scale increment $scaleincrementtomeetmeanharmony")
        display(additionalcapital)
    end
    # now find which of the previous years will supply this
    yearsPrevious = [y for y = 1:destyear-1]
    bestyear = 0
    bestgain = 0
    bestscenario = s
    for y in yearsPrevious
        originalcapital = additionalcapital
        if inversedepreciateinvestments
            originalcapital =
                inversedepreciate(additionalcapital, destyear - 1 - y, s.prob.dep)
        end
        gain, newscenario = gainfromInvesting(s, y, originalcapital)# negative
        posflags = (s.netoutputs[y][1:n] .> zeros(n))

        if (gain) > bestgain && *(posflags...)
            bestyear = y
            bestgain = gain
            bestscenario = newscenario
        end
    end
    if bestyear == 0
        # we have not found a possible transfer of resources
        # print warning and terminate
        global doagain = false
        println(" we have not found a possible transfer of resources, terminating")
        return s
    end
    if verbose
        println(
            " \nto achieve gain ",
            bestgain,
            " take from year ",
            bestyear,
            " lowest harmony year was=",
            destyear,
            "\n harmonies =",
            s.h,
            "\n capital=",
            capitalbyyear(s),
            "\n new harmonies =",
            bestscenario.h,
            "\n new capital ",
            capitalbyyear(bestscenario),
        )

    end
    return bestscenario
end
function Estimate_how_much__production_to_be_scaled_up(s::Scenario, this_year)
    target = harmonyInverse(s.meanh)
    current = harmonyInverse(s.h[this_year])
    diff = target - current
    if diff < 0
        println("\nrequest for negative upscale $diff, target=$target, current=$current ")
        println(
            "\nharmony for year $this_year was ",
            s.h[this_year, :],
            " mean was ",
            s.meanh,
        )
        throw(error("request for negative upscale"))
    end
    fractiontomove = diff * epsilon
    if verbose
        println("\nupscale ", this_year, " by ", fractiontomove)
    end
    return fractiontomove

end
function readInProblem(flowname, capname, depname, labtagsname)
    headers, flows, caps, dep, labtarg =
        readinspreadsheets(flowname, capname, depname, labtagsname)
    # Generate key matrices
    A = generateAMatrix(flows)
    V = compute_Labour_Value_From_Augmented_A_Matrix(A)

    C = generateCMatrix(flows, caps)
    D = generateDMatrix(C, dep)

    Horizons = 1 ./ dep

    AD = A + D
    IminusA = I - A

    # process the targets and labour supply
    labouravailable = (extractlaboursupply(labtarg, headers))
    g = extracttargets(labtarg, headers)
    TheLastYear, gcol = size(g)
    caprows, capcols = size(caps)
    if verbose
        println("A")
        display(A)
        println()
        println("D")
        display(D)
        println()
        println("inv(I-A)")
        display(inv(I - A))
        println("labour")
        display(labouravailable)
        println()
        display(headers)
        println()

        println("net consumer demand goals")
        display(g)
        println()
    end
    return PlanProblem(
        headers,
        flows,
        caps,
        dep,
        labtarg,
        A,
        V,
        C,
        D,
        AD,
        IminusA,
        inv(IminusA),
        labouravailable,
        g,
        TheLastYear,
        caprows,
        capcols,
    )
end
function Compute_the_degree_of_goal_fulfillment_that_is_in_principle_possible_goal_fullfilment_ratio_vector(
    s::Scenario,
)
    s.goal_fullfilment_ratio_vector =
        Compute_the_degree_of_goal_fulfillment_that_is_in_principle_possible_goal_fullfilment_ratio_vector(
            s.si,
            s.prob.labouravailable,
            s.O,
            s.prob.A,
            s.prob.C,
        )
    # println("goal fullfilment ratio ",s.goal_fullfilment_ratio_vector)
end

#------------------------------------ end of pure functions
function solvePlanProblem(fln, cpn, dpn, ltn)
    #------------------------------------ some   vars for main function
    iter = 0

    #- ---------------------------------- Executable statements below

    #=
    	read in the files from cvs format
    	=#
    problem = readInProblem(fln, cpn, dpn, ltn)
    #=

    	1 For the last year of the plan set a net output target such that \label{item:settarget}

    	    a Gross output is such as to ensure full employment of the workforce
    	    b Sufficient investment is being carried out to compensate for depreciation during the year.

    	 This net output target must be a scaled version of the originally specified
    	 target. If the original target would have caused unemployment it will be up-scaled
    	 otherwise down-scaled.=#

    lastinvest, lastgoal = For_the_last_year_of_the_plan_return_a_net_output_target(
        problem.TheLastYear,
        problem.A,
        problem.D,
        problem.g,
        problem.labouravailable,
    )
    if verbose
        display(lastinvest)
        display(lastgoal)
    end
    problem.g[problem.TheLastYear, :] = lastgoal

    #= 1.1 Also compute the gross outputs required each year by the net output
    	       goals =#
    Otmp =
        ([grossOutputForDemandf(problem.g[i, :], problem.A)' for i = 1:problem.TheLastYear])

    # make sure it is a matrix
    Otmp = vcat(Otmp...)

    #=
    	2 Assign to each year's fixed capital stock the starting stock depreciated by
    	 the appropriate depreciation rate on that type of fixed capital stock.
    	 This is a rank 3 tensor indexed by year, source industry, destination industry.
    	 Create an initially zeroed raank 3 tensor of investments and add to it a preliminary investment schedule
    	 =#
    sitmp =
        Assign_to_each_year_capital_stock(problem.caps, problem.dep, problem.TheLastYear)
    investmentstmp = zeros(Float64, problem.TheLastYear, problem.caprows, problem.capcols)
    preassignedcapital = initialinvestmentlevel * (problem.caps .* problem.dep)
    for y = 1:problem.TheLastYear-1
        setup_preliminary_investment_schedule(
            y,
            investmentstmp,
            preassignedcapital,
            problem.dep,
            sitmp,
        )
    end


    baseScenario = Scenario(
        problem,
        Otmp,
        sitmp,
        investmentstmp,
        zeros(Float64, problem.TheLastYear), #goal fullfillment
        zeros(Float64, problem.TheLastYear), #h
        0,# meanh
        0,# stdevh
        zeros(Float64, problem.TheLastYear, problem.capcols),
        zeros(Float64, problem.TheLastYear, problem.capcols),
        deepcopy(problem.g),
    )
    #=
              3 Compute the degree of goal fulfillment that is in principle possible
              given these fixed capital stocks and the available workforce. From this compute
              the Harmonies.
              =#
    #= 3.1 Also compute the gross outputs required each year by the net output
           goals this time including investment in the goals=#
    update_outputs(baseScenario)
    #=
     4 Compute the mean harmony, and standard deviation of the harmonies over the whole time period.
      =#

    computeHarmonies(baseScenario)
    iter = 1
    while doagain
        for i = 2:problem.TheLastYear


            #=
             5 If the coefficient of harmony variation has fallen below some threshold, terminate the plan prepartion, otherwise continue to   step 6.
             =#
            if ((baseScenario.stdh / abs(baseScenario.meanh)) < mincoeff) ||
               (iter > maxiter)
                println(
                    " iterations,Harmony coefficient of variation \n",
                    iter,
                    ",",
                    baseScenario.stdh / baseScenario.meanh,
                )
                global doagain = false
                break
            end
            iter = iter + 1
            if doagain
                #=
                6 Select a year with  low harmony.
                 	 =#

                if baseScenario.h[i] < baseScenario.meanh
                    # println("$iter,$i, ",baseScenario.h[i] ," ", baseScenario.meanh)
                    lowyear = i  #Select_the_year_with_the_lowest_harmony(baseScenario.h)
                    #=
                    7 Estimate by how much the  production  of this year would have to be scaled up in order to be at the mean harmony.
                    This is returned as a ratio.
                    To do this we first work out the level of net production needed to meet the mean harmony. This involves inverting the harmony function.
                     	=#
                    upscale =
                        Estimate_how_much__production_to_be_scaled_up(baseScenario, lowyear)
                    #=
                    8 Attempt to scale up some fraction $\epsilon$ of the way towards the mean harmony by \label{item:realloc}
                    scheduling investment in  at least one previous year. This will obviously alter both the consumption in previous years
                    and the fixed capital stocks in one or more years.
                             	 =#

                    baseScenario = Attempt_to_scale_up(baseScenario, lowyear, upscale)
                end
            end
            #=
               9 Go back to step 5
             =#
        end
        if !doagain
            for yr = 1:problem.TheLastYear-depreciationhorizon
                display_result(
                    yr,
                    baseScenario.investmentsByTypeAndYear,
                    problem.C,
                    baseScenario.si,
                    baseScenario.h,
                    problem.headers,
                    baseScenario.O[yr, :] .* baseScenario.goal_fullfilment_ratio_vector[yr],
                    problem.A,
                    problem.labouravailable[yr],
                    problem.g,
                    baseScenario,
                )
            end
        end
    end
end

solvePlanProblem("jeuflows.csv", "jeucap.csv", "jeudep.csv", "jeulabtargs.csv")

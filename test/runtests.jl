using Test
using CoNCMOR

module mpart3
using StaticArrays
using SparseArrays
using Test
# using PlotlyJS
using FinEtools
using CoNCMOR: CoNCData, nclusters, nfuncspercluster
using CoNCMOR: nbasisfunctions, transfmatrix, LegendreBasis
colors = ["rgb(164, 194, 244)", "rgb(194, 194, 144)", "rgb(194, 144, 244)", "rgb(164, 244, 144)", "rgb(164, 194, 244)", "rgb(255, 217, 102)", "rgb(234, 153, 153)", "rgb(142, 124, 195)"]
function test()
    points = StaticArrays.SArray{Tuple{2},Float64,1,2}[[0.6501312860170676, 0.889978628556229], [0.19162369745368268, 0.5000498373796463], [0.9513777178320113, 0.9313997945413548], [0.6715708752191747, 0.0031159439825900748], [0.908983887358429, 0.33697103964932174], [0.8630730459099334, 0.8172618688642774], [0.21191150611345178, 0.6279042354095665], [0.0017150255179243512, 0.7401108604352873], [0.045937411259845184, 0.3867144891772889], [0.07794203801162314, 0.45466525770422384], [0.11103905573968631, 0.7213589211503646], [0.522501483190954, 0.9505399568028736], [0.14011031960881337, 0.5976293299997106], [0.15824295891177775, 0.4538838028009915], [0.22532364404660643, 0.7490903677067646], [0.16495650937348727, 0.4095197672868107], [0.2929326909660668, 0.2627344309574493], [0.7749630590488474, 0.010281406373283675], [0.23137166232629092, 0.1872554576301051], [0.48014323786400226, 0.8846143614245845], [0.26528544579757507, 
    0.6924068278632542], [0.2839884809072999, 0.09439141875666568], [0.942287809602099, 0.4005321413670275], [0.026409937646382442, 0.6845508696068237], [0.7752785557256849, 0.3285662200502397], [0.09380966727246509, 0.15933031973012612], [0.2813301059379758, 0.15216808040457552], [0.4780694864591679, 0.6063568972874278], [0.6982811071004462, 0.5308171330635794], [0.9609668367122937, 0.09387988328210572]]  
    X = fill(0.0, length(points), 2)
    for i in 1:length(points)
        X[i, :] .= points[i]
    end        
    ppartitioning = pointpartitioning(X, 2)
    partitionnumbers = unique(ppartitioning)

    # data = PlotlyBase.AbstractTrace[]
    # for gp in partitionnumbers
    #     trace1 = scatter(; 
    #         x=[points[i][1] for i in 1:length(points) if ppartitioning[i] == gp], 
    #         y=[points[i][2] for i in 1:length(points) if ppartitioning[i] == gp] ,
    #       mode="markers",
    #       marker=attr(color=colors[gp], size=12,
    #           line=attr(color="white", width=0.5))
    #       )
    #     push!(data, trace1)
    # end
    # layout = Layout(;title="Point partitioning",
    #     xaxis=attr(title="x", zeroline=false),
    #     yaxis=attr(title="y", zeroline=false))

    # pl = plot(data, layout)
    # display(pl)

    mor = CoNCData(X, ppartitioning)
    @test nclusters(mor) == length(partitionnumbers)
    geom = NodalField(X) 
    u = NodalField(X) 
    numberdofs!(u)
    Phi  = transfmatrix(mor, LegendreBasis, 1, u)
    @test nfuncspercluster(mor) == 1
    @test nbasisfunctions(mor) == 2
    
    I, J, V = findnz(Phi)
    @test (I, J, V) == ([1, 5, 7, 9, 11, 23, 35, 45, 49, 55, 57, 59, 2, 6, 8, 10, 12, 24, 36, 46, 50, 56, 58, 60, 3, 13, 15, 17, 19, 21, 25, 27, 29, 31, 33, 37, 39, 41, 43, 47, 51, 53, 4, 14, 16, 18, 20, 22, 26, 28, 30, 32, 34, 38, 40, 42, 44, 48, 52, 54], [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 
    4, 4, 4, 4, 4, 4, 4], [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0])
    nothing
end
end
using .mpart3
mpart3.test()



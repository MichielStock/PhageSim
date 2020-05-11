#Initiëren
start = time()
nsims = 1000
eindbact = zeros(Int64,3,1)
tmpbact = zeros(Int64,3,1)
eindfag = zeros(Int64,3,1)
tmpfag = zeros(Int64,3,1)

#Simulaties nsims keer doorlopen en de grootte van de bacteriële populatie
#per species op de laatste tijdstap teruggeven
for i = 1:nsims
    include("scripts/shortbachelorprojectsim.jl")
    tmpbact[1:3] = last(results)[1]
    global eindbact = hcat(eindbact,tmpbact)
end

eindbact = eindbact[:,2:end]

#Kijken hoeveel populaties in evenwicht zijn
nevenwicht = sum(eindbact .!=0,dims=1)

nul = 0
een = 0
twee = 0
drie = 0
for j = 1:length(nevenwicht)
    if nevenwicht[j] == 0
        global nul = nul+1
    elseif nevenwicht[j] == 1
        global een = een+1
    elseif nevenwicht[j] == 2
        global twee = twee+1
    else
        global drie = drie+1
    end
end

resultaat = hcat(drie,twee,een,nul)
resultaatproc = resultaat./sum(resultaat)

elapsed = time() - start

# setup CMSSW to run this code, eg:
#  cmsrel CMSSW_12_0_2 && cd CMSSW_12_0_2/src && cmsenv

import ROOT as r
from glob import glob
from array import array
# need fwlite to resolve the references
from DataFormats.FWLite import Events, Handle

out = r.TFile('qg.root', 'recreate')

tjet = r.TTree('jets', 'jets')
quark = array('i', [0]); tjet.Branch('quark', quark, 'quark/I')
nconstituents = array('i', [0]); tjet.Branch('nconstituents', nconstituents, 'nconstituents/I')
eventn = array('i', [0]); tjet.Branch('eventn', eventn, 'eventn/I')
pt = array('f', [0]); tjet.Branch('pt', pt, 'pt/F')
eta = array('f', [0]); tjet.Branch('eta', eta, 'eta/F')
phi = array('f', [0]); tjet.Branch('phi', phi, 'phi/F')
mass = array('f', [0]); tjet.Branch('mass', mass, 'mass/F')
area = array('f', [0]); tjet.Branch('area', area, 'area/F')
ptD = array('f', [0]); tjet.Branch('ptD', ptD, 'ptD/F')

tevt = r.TTree('event', 'event')
npartons = array('i', [0]); tevt.Branch('npartons', npartons, 'npartons/I')
nrecojets = array('i', [0]); tevt.Branch('nrecojets', nrecojets, 'nrecojets/I')

#t = r.TChain('Events')
d='/eos/cms/store/mc/RunIISummer20UL18MiniAODv2/QCD_HT100to200_TuneCP5_PSWeights_13TeV-madgraph-pythia8/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/2520000'
#_=[t.Add(f) for f in glob(d+'/*')[:2]]
f = r.TFile(glob(d+'/*')[0])

ntot, ngen, ngood = 0,0,0
for e in f.Events:
    ntot += 1
    # Find partons
    ppdg = []; pp4 = []
    for p in f.Events.recoGenParticles_prunedGenParticles__PAT.product():
        if p.status() != 23: continue
        ppdg.append(p.pdgId())
        pp4.append(p.p4())
    # lets avoid events with heavy flavour
    if 4 in ppdg or 5 in ppdg or -4 in ppdg or -5 in ppdg: continue
    njets = len(ppdg)
    # have the genjet topology match the partons
    if f.Events.recoGenJets_slimmedGenJets__PAT.size() != njets: continue
    # find genjets associated to partons
    match = [False]*njets
    g2p = [-1]*njets
    for igj, gjet in enumerate(f.Events.recoGenJets_slimmedGenJets__PAT.product()):
        best = -1; best_dr = 9999.
        jet4 = gjet.p4()
        for ipart, p4 in enumerate(pp4):
            dr = r.Math.VectorUtil.DeltaR(jet4,p4)
            if dr < best_dr:
                best_dr = dr
                best = ipart
        if best_dr < 0.4:
            match[best] = True
            g2p[igj] = best
    if sum(match) == njets:
        ngen += 1
    else:
        continue    
    # check reco jets
    nreco = f.Events.patJets_slimmedJets__PAT.size()
    g2j = [-1]*njets
    match = [False]*njets
    for ij, jet in enumerate(f.Events.patJets_slimmedJets__PAT.product()):
        best = -1; best_dr = 9999.
        jet4 = jet.p4()
        for igj, gjet in enumerate(f.Events.recoGenJets_slimmedGenJets__PAT.product()):
            dr = r.Math.VectorUtil.DeltaR(jet4,gjet.p4())
            if dr < best_dr:
                best_dr = dr
                best = igj
        if best_dr < 0.4:
            match[best] = True
            g2j[best] = ij
    # events with one reco jet matching to one gen jet matching to one parton
    # allow more reco jets since we can have pileup
    if sum(match) == njets:
        ngen += 1
    else:
        continue    
    ngood += 1
    eventn[0] = tevt.GetEntries()
    npartons[0] = njets
    nrecojets[0] = f.Events.patJets_slimmedJets__PAT.size()
    _ = tevt.Fill()
    for igj,ij in enumerate(g2j):
        jet = f.Events.patJets_slimmedJets__PAT.at(ij)
        quark[0] = 0 if f.Events.recoGenParticles_prunedGenParticles__PAT.at(g2p[igj]).pdgId() == 22 else 1
        pt[0] = jet.pt()
        eta[0] = jet.eta()
        phi[0] = jet.phi()
        mass[0] = jet.mass()
        area[0] = jet.jetArea()
        sum_pt = 0.
        sum_weight = 0.
        for daught in jet.getJetConstituents():
            sum_pt += daught.pt()
            sum_weight += daught.pt()*daught.pt()
        ptD[0] = r.TMath.Sqrt(sum_weight) / sum_pt;
        _ = tjet.Fill()

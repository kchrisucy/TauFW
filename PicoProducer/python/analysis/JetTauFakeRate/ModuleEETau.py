# Author: Konstantinos Christoforou (Jan 2022)
# Description: Simple module to pre-select eetau events

import sys
import numpy as np

from ROOT import TFile, TTree, TH1D
#from TauFW.PicoProducer.analysis.ModuleTauTriad import *
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from TauPOG.TauIDSFs.TauIDSFTool import TauIDSFTool, TauESTool, campaigns

from TauFW.PicoProducer import datadir
#from TauFW.PicoProducer.analysis.utils import *
from TauFW.PicoProducer.corrections.MuonSFs import *
from TauFW.PicoProducer.corrections.BTagTool import BTagWeightTool, BTagWPs
# from TauFW.PicoProducer.analysis.utils import ensurebranches, redirectbranch, deltaPhi, getmet, getmetfilters, correctmet, getlepvetoes
from TauFW.PicoProducer.analysis.utils import deltaPhi, getmetfilters
from TauFW.PicoProducer.corrections.PileupTool import *
#from TauFW.PicoProducer.corrections.MuonSFs import *
from TauFW.PicoProducer.corrections.ElectronSFs import *
from TauFW.PicoProducer.corrections.TrigObjMatcher import TrigObjMatcher
from TauPOG.TauIDSFs.TauIDSFTool import TauIDSFTool, TauESTool, campaigns

class ModuleEETau(Module):
  
  def __init__(self,fname,**kwargs):
    self.outfile = TFile(fname,'RECREATE')
    self.channel    = kwargs.get('channel',  'none'         ) # channel name
    self.dtype      = kwargs.get('dtype',    'data'         )
    self.ismc       = self.dtype=='mc'
    self.isdata     = self.dtype=='data' or self.dtype=='embed'
    self.isembed    = self.dtype=='embed'
    self.era        = kwargs.get('era',     2017           ) # integer, e.g. 2017, 2018, 2018UL
    #self.dojec      = kwargs.get('jec',      True           ) and self.ismc #and self.year==2016 #False
    #self.useT1      = kwargs.get('useT1',    False          ) # MET T1
    self.verbosity  = kwargs.get('verb',     0              ) # verbosity
    
    tauSFVersion  = { 2016: '2016Legacy', 2017: '2017ReReco', 2018: '2018ReReco' }
    if "2016" in self.era : 
      self.yearForTauSF = 2016
    elif "2017" in self.era : 
      self.yearForTauSF = 2017
    elif "2018" in self.era :
      self.yearForTauSF = 2018
    else :
      "Sanity check error, no 2016/2017/2018 in self.era"

    jsonfile       = os.path.join(datadir,"trigger/tau_triggers_%d.json"%(self.yearForTauSF))
    self.trigger   = TrigObjMatcher(jsonfile,trigger='SingleElectron',isdata=self.isdata)

    ## CORRECTIONS
    self.metfilt      = getmetfilters(self.era,self.isdata,verb=self.verbosity)
    if self.ismc:
      self.puTool     = PileupWeightTool(era=self.era,sample=fname,verb=self.verbosity)
      self.eleSFs     = ElectronSFs(era=self.era) # electron id/iso/trigger SFs
      self.tesTool    = TauESTool(tauSFVersion[self.yearForTauSF])  # real tau energy scale corrections
      self.tauSFsT_dm = TauIDSFTool(tauSFVersion[self.yearForTauSF],'DeepTau2017v2p1VSjet','Tight', dm=True)
      self.tauSFsT    = TauIDSFTool(tauSFVersion[self.yearForTauSF],'DeepTau2017v2p1VSjet','Tight')
      self.tauSFsM    = TauIDSFTool(tauSFVersion[self.yearForTauSF],'DeepTau2017v2p1VSjet','Medium')
      self.tauSFsVVVL = TauIDSFTool(tauSFVersion[self.yearForTauSF],'DeepTau2017v2p1VSjet','VVVLoose')      
      self.etfSFs     = TauIDSFTool(tauSFVersion[self.yearForTauSF],'DeepTau2017v2p1VSe',  'Medium')
      self.mtfSFs     = TauIDSFTool(tauSFVersion[self.yearForTauSF],'DeepTau2017v2p1VSmu', 'Tight')
    ## after opening the root files for muonSF and tauSF, something is getting confused, so you should reopen
    ## your output file in order to store your trees/branches
    self.outfile.cd()

    
    #print "Hi you are in mumutau after muSF"

    self.jetCutPt   = 30
    self.bjetCutEta = 2.4 if self.era==2016 else 2.5
    self.deepjet_wp = BTagWPs('DeepJet',era=self.era)

    # TRIGGER THRESHOLDS
    self.triggerPtThreshold = 0.0
    if self.era=="2016" or self.era=="UL2016":
      self.triggerPtThreshold = 30.0
    elif self.era=="2017" or self.era=="UL2017" or self.era=="2018" or self.era=="UL2018":
      self.triggerPtThreshold =35.0

  def beginJob(self):
    """Prepare output analysis tree and cutflow histogram."""
    
    # CUTFLOW HISTOGRAM
    self.cutflow  = TH1D('cutflow','cutflow',25,0,25)
    self.cut_none     = 0
    self.cut_trig     = 1
    self.cut_electron = 2
    self.cut_tau      = 3
    self.cut_pair     = 4
    self.cutflow.GetXaxis().SetBinLabel(1+self.cut_none,     "no cut"      )
    self.cutflow.GetXaxis().SetBinLabel(1+self.cut_trig,     "trigger"     )
    self.cutflow.GetXaxis().SetBinLabel(1+self.cut_electron, "electron"    )
    self.cutflow.GetXaxis().SetBinLabel(1+self.cut_tau,      "tau"         )
    self.cutflow.GetXaxis().SetBinLabel(1+self.cut_pair,     "pair"        )


    # TREE
    self.tree   = TTree('tree','tree')
    self.evt           = np.zeros(1,dtype='i')
    self.lumi          = np.zeros(1,dtype='i')
    self.metfilter     = np.zeros(1,dtype='?')
    self.genweight     = np.zeros(1,dtype='f')
    self.puweight      = np.zeros(1,dtype='f')
    self.trigweight    = np.zeros(1,dtype='f')
    self.idisoweight_1 = np.zeros(1,dtype='f')
    self.idisoweight_2 = np.zeros(1,dtype='f')
    self.idweightTdm_tau  = np.zeros(1,dtype='f')
    self.idweightT_tau    = np.zeros(1,dtype='f')
    self.idweightM_tau    = np.zeros(1,dtype='f')
    self.idweightVVVL_tau = np.zeros(1,dtype='f')
    self.ltfweight_tau    = np.zeros(1,dtype='f')
    self.ltfweight_tau    = np.zeros(1,dtype='f')
    self.pt_e0   = np.zeros(1,dtype='f')
    self.eta_e0  = np.zeros(1,dtype='f')
    self.q_e0    = np.zeros(1,dtype='i')
    self.id_e0   = np.zeros(1,dtype='?')
    self.iso_e0  = np.zeros(1,dtype='f')
    self.pt_e1   = np.zeros(1,dtype='f')
    self.eta_e1  = np.zeros(1,dtype='f')
    self.q_e1    = np.zeros(1,dtype='i')
    self.id_e1   = np.zeros(1,dtype='?')
    self.iso_e1  = np.zeros(1,dtype='f')
    self.pt_tau   = np.zeros(1,dtype='f')
    self.eta_tau  = np.zeros(1,dtype='f')
    self.q_tau    = np.zeros(1,dtype='i')
    self.id_tau   = np.zeros(1,dtype='?')
    self.iso_tau  = np.zeros(1,dtype='f')
    ## Jet to tau FR
    self.IsOnZ              = np.zeros(1,dtype='?')
    self.TauIsGenuine       = np.zeros(1,dtype='?')
    self.TauPt              = np.zeros(1,dtype='f')
    self.TauEta             = np.zeros(1,dtype='f')
    self.TauDM              = np.zeros(1,dtype='f')
    self.JetN               = np.zeros(1,dtype='i')
    self.BJetN              = np.zeros(1,dtype='i')
    self.HT                 = np.zeros(1,dtype='f')
    self.LeptonOnePt        = np.zeros(1,dtype='f')
    self.LeptonTwoPt        = np.zeros(1,dtype='f')
    self.DileptonPt         = np.zeros(1,dtype='i')
    self.DileptonMass       = np.zeros(1,dtype='f')
    self.DileptonDeltaEta   = np.zeros(1,dtype='f')
    self.DileptonDeltaPhi   = np.zeros(1,dtype='f')
    self.DileptonDeltaR     = np.zeros(1,dtype='f')

    self.tree.Branch('evt',          self.evt,           'evt/I'          )          
    self.tree.Branch('lumi',         self.lumi,          'lumi/F'         )
    self.tree.Branch('metfilter',    self.metfilter,     'metfilter/O'    )              
    self.tree.Branch('genweight',    self.genweight,     'genweight/F'    )          
    self.tree.Branch('puweight',     self.puweight,      'puweight/F'     )          
    self.tree.Branch('trigweight',   self.trigweight,    'trigweight/F'   )
    self.tree.Branch('idisoweight_1',self.idisoweight_1, 'idisoweight_1/F')
    self.tree.Branch('idisoweight_2',self.idisoweight_2, 'idisoweight_2/F')
    self.tree.Branch('idweightTdm_tau',  self.idweightTdm_tau,  'idweightTdm_tau/F'  )
    self.tree.Branch('idweightT_tau',    self.idweightT_tau,    'idweightT_tau/F'    )
    self.tree.Branch('idweightM_tau',    self.idweightM_tau,    'idweightM_tau/F'    )
    self.tree.Branch('idweightVVVL_tau', self.idweightVVVL_tau, 'idweightVVVL_tau/F' )
    self.tree.Branch('ltfweight_tau',    self.ltfweight_tau,    'ltfweight_tau/F'    ) 
    self.tree.Branch('ltfweight_tau',    self.ltfweight_tau,    'ltfweight_tau/F'    )
    self.tree.Branch('pt_e0',   self.pt_e0,  'pt_e0/F' )
    self.tree.Branch('eta_e0',  self.eta_e0, 'eta_e0/F')
    self.tree.Branch('q_e0',    self.q_e0,   'q_e0/I'  )
    self.tree.Branch('id_e0',   self.id_e0,  'id_e0/O' )
    self.tree.Branch('iso_e0',  self.iso_e0, 'iso_e0/F')
    self.tree.Branch('pt_e1',   self.pt_e1,  'pt_e1/F' )
    self.tree.Branch('eta_e1',  self.eta_e1, 'eta_e1/F')
    self.tree.Branch('q_e1',    self.q_e1,   'q_e1/I'  )
    self.tree.Branch('id_e1',   self.id_e1,  'id_e1/I' )
    self.tree.Branch('iso_e1',  self.iso_e1, 'iso_e1/F')
    self.tree.Branch('pt_tau',   self.pt_tau,  'pt_tau/F' )
    self.tree.Branch('eta_tau',  self.eta_tau, 'eta_tau/F')
    self.tree.Branch('q_tau',    self.q_tau,   'q_tau/I'  )
    self.tree.Branch('id_tau',   self.id_tau,  'id_tau/I' )
    self.tree.Branch('iso_tau',  self.iso_tau, 'iso_tau/F')
    ## Jet to tau FR
    self.tree.Branch("IsOnZ"           ,  self.IsOnZ           , "IsOnZ/O"            )           
    self.tree.Branch("TauIsGenuine"    ,  self.TauIsGenuine    , "TauIsGenuine/O"     )           
    self.tree.Branch("TauPt"           ,  self.TauPt           , "TauPt/F"            )           
    self.tree.Branch("TauEta"          ,  self.TauEta          , "TauEta/F"           ) 
    self.tree.Branch("TauDM"           ,  self.TauDM           , "TauDM/I"            ) 
    self.tree.Branch("JetN"            ,  self.JetN            , "JetN/I"             ) 
    self.tree.Branch("BJetN"           ,  self.BJetN           , "BJetN/I"            ) 
    self.tree.Branch("HT"              ,  self.HT              , "HT/F"               ) 
    self.tree.Branch("LeptonOnePt"     ,  self.LeptonOnePt     , "LeptonOnePt/F"      ) 
    self.tree.Branch("LeptonTwoPt"     ,  self.LeptonTwoPt     , "LeptonTwoPt/F"      ) 
    self.tree.Branch("DileptonPt"      ,  self.DileptonPt      , "DileptonPt/F"       ) 
    self.tree.Branch("DileptonMass"    ,  self.DileptonMass    , "DileptonMass/F"     ) 
    self.tree.Branch("DileptonDeltaEta",  self.DileptonDeltaEta, "DileptonDeltaEta/F" )
    self.tree.Branch("DileptonDeltaPhi",  self.DileptonDeltaPhi, "DileptonDeltaPhi/F" ) 
    self.tree.Branch("DileptonDeltaR"  ,  self.DileptonDeltaR  , "DileptonDeltaR/F"   )

                                            
  def endJob(self):                         
    """Wrap up after running on all events and files"""
    self.outfile.Write()
    self.outfile.Close()
  
  def analyze(self, event):
    """Process event, return True (pass, go to next module) or False (fail, go to next event)."""
    
    # EVENT
    #print "New event"
    #print event.event
    #eventNum = event.event & 0xffffffffffffffff
    #self.evt[0]             = eventNum
    #print type(event.event)
    #print event.event & 0xffffffffffffffff #, type(event.event & 0xffffffffffffffff)
    self.evt[0]             = event.event & 0xffffffffffffffff 
    # NO CUT
    self.cutflow.Fill(self.cut_none)
    
    # TRIGGER
    if self.era=="2016" or self.era=="UL2016":
      if not (event.HLT_Ele27_WPTight_Gsf) : return False
    elif self.era=="2017" or self.era=="UL2017":
      if not (event.HLT_Ele32_WPTight_Gsf_L1DoubleEG and event.TrigObj_filterBits == 1024) : return False # 32_L1DoubleEG_AND_L1SingleEGOr filter to extract ele27 from ele35
    elif self.era=="2018" or self.era=="UL2018":
      if not event.HLT_Ele32_WPTight_Gsf : return False
    self.cutflow.Fill(self.cut_trig)
    
    # SELECT ELECTRONS
    electrons = [ ]
    for electron in Collection(event,'Electron'):
      if electron.pt<10: continue
      if abs(electron.eta)>2.1: continue
      if abs(electron.dz)>0.2: continue
      if abs(electron.dxy)>0.045: continue
      if electron.lostHits>1: continue
      if not (electron.mvaFall17V2Iso_WP90) : continue
      if electron.pfRelIso03_all > 0.50: continue
      electrons.append(electron)
    if len(electrons)!=2: return False
    
    passTriggerMatching = False
    for selectedElectron in electrons:
      if electron.pt < self.triggerPtThreshold : continue
      if not self.trigger.match(event,electron): continue
      passTriggerMatching = True
    #kc
    if not passTriggerMatching : return False
    ###
    self.cutflow.Fill(self.cut_electron)
    
    # SELECT TAU
    taus = [ ]
    for tau in Collection(event,'Tau'):
      if tau.pt<20: continue
      if abs(tau.eta)>2.4: continue
      if abs(tau.dz)>0.2: continue
      if tau.decayMode not in [0,1,10,11]: continue
      if abs(tau.charge)!=1: continue
      if tau.idDeepTau2017v2p1VSe<16: continue # medium Vse
      if tau.idDeepTau2017v2p1VSmu<8: continue # tight Vsmu 
      if tau.idDeepTau2017v2p1VSjet<1: continue # start with VVVL versusJets
      taus.append(tau)
    if len(taus)!=1: return False
    self.cutflow.Fill(self.cut_tau)

    # Leptons
    electron0 = electrons[0]
    electron1 = electrons[1]
    tau       = max(taus,key=lambda p: p.pt)

    # Leptons dR
    #muon = max(muons,key=lambda p: p.pt)
    if electron0.DeltaR(electron1)<0.4: return False
    if electron0.DeltaR(tau)<0.4: return False
    if electron1.DeltaR(tau)<0.4: return False
    self.cutflow.Fill(self.cut_pair)
    
    #KC checks
    #for i, muon1 in enumerate(muons,1):
    #  print "The %s muon has pT = %s" %(i, muon1.pt)
    #  
    #print "OK, let's see this"
    #print "muon = max(muons,key=lambda p: p.pt), pT = %s" %(muon.pt)
    #print "muon0 = muons[0], pT = %s" %(muon0.pt)
    #print "muon1 = muons[1], pT = %s" %(muon1.pt)

    # VETO for extraElectrons
    extraelectron_veto = False
    looseElectrons = [ ]    
    for electron in Collection(event,'Electron'):
      if electron.pt<10: continue
      if abs(electron.eta)>2.4: continue
      if abs(electron.dz)>0.2: continue
      if abs(electron.dxy)>0.045: continue
      if electron.lostHits>1: continue
      if electron.pfRelIso03_all>0.3: continue
      if any(electron.DeltaR(tau)<0.4 for tau in taus): continue
      if electron.mvaFall17V2noIso_WPL and all(e._index!=electron._index for e in electrons):
        looseElectrons.append(electron)
        extraelectron_veto = True
        
    if extraelectron_veto : return False
    
    # VETO MUONS
    muon_veto = False

    looseMuons = [ ]
    for muon in Collection(event,'Muon'):
      if muon.pt<10: continue
      if abs(muon.eta)>2.5: continue
      if abs(muon.dz)>0.2: continue
      if abs(muon.dxy)>0.045: continue
      if muon.pfRelIso03_all>0.3: continue
      if any(muon.DeltaR(tau)<0.4 for tau in taus): continue
      if muon.convVeto==1 and muon.lostHits<=1 and muon.mvaFall17V2noIso_WPL:
        looseMuons.append(muon)
        muon_veto = True
        
    if muon_veto : return False
    
    # SELECT JETS
    #jets, bjets, met = fillJetBranches(self,event,tau)
    jets, bjets = fillJetBranches(self,event,tau)
    is0b = False
    if(len(bjets) == 0) : is0b = True
    MET = event.MET_pt    
    MET_phi = event.MET_phi

    # WEIGHTS AND CORRECTIONS
    # EVENT WEIGHTS
    self.lumi[0]            = event.luminosityBlock
    self.metfilter[0]       = self.metfilt(event)
    
    if self.ismc : 
      self.genweight[0]     = event.genWeight
      self.puweight[0]      = self.puTool.getWeight(event.Pileup_nTrueInt)
      # ELECTRON WEIGHTS
      self.trigweight[0]    = self.eleSFs.getTriggerSF(electron0.pt,electron0.eta) # assume leading electron was triggered on
      self.idisoweight_1[0] = self.eleSFs.getIdIsoSF(electron0.pt,electron0.eta)
      self.idisoweight_2[0] = self.eleSFs.getIdIsoSF(electron1.pt,electron1.eta)
      # TAU WEIGHTS
      self.idweightTdm_tau[0]  = 1.0
      self.idweightT_tau[0]    = 1.0
      self.idweightM_tau[0]    = 1.0
      self.idweightVVVL_tau[0] = 1.0
      self.ltfweight_tau[0]    = 1.0
      self.ltfweight_tau[0]    = 1.0
      
      if tau.genPartFlav==5: # real tau
        self.idweightTdm_tau[0]  = self.tauSFsT_dm.getSFvsDM(tau.pt,tau.decayMode)
        self.idweightT_tau[0]    = self.tauSFsT.getSFvsPT(tau.pt)
        self.idweightM_tau[0]    = self.tauSFsM.getSFvsPT(tau.pt) 
        self.idweightVVVL_tau[0] = self.tauSFsVVVL.getSFvsPT(tau.pt)
      elif tau.genPartFlav in [1,3]: # e -> tau fake                       
        self.ltfweight_tau[0]    = self.etfSFs.getSFvsEta(tau.eta,tau.genPartFlav)
      elif tau.genPartFlav in [2,4]: # mu -> tau fake                             
        self.ltfweight_tau[0]    = self.mtfSFs.getSFvsEta(tau.eta,tau.genPartFlav)


    # SAVE CONTROL VARIABLES
    self.pt_e0[0]   = electron0.pt
    self.eta_e0[0]  = electron0.eta
    self.q_e0[0]    = electron0.charge
    self.id_e0[0]   = electron0.cutBased
    self.iso_e0[0]  = electron0.pfRelIso03_all
    self.pt_e1[0]   = electron1.pt
    self.eta_e1[0]  = electron1.eta
    self.q_e1[0]    = electron1.charge
    self.id_e1[0]   = electron1.cutBased
    self.iso_e1[0]  = electron1.pfRelIso03_all
    self.pt_tau[0]  = tau.pt
    self.eta_tau[0] = tau.eta
    self.q_tau[0]   = tau.charge
    self.id_tau[0]  = tau.idDeepTau2017v2p1VSjet
    self.iso_tau[0] = tau.rawIso

    #KC checks
    #print "Hi there, you are in ModuleMuMuTau, under #KC checks-----"
    #for i, jet in enumerate(jets,1):
    #  print "The %s jet has pT = %s" %(i, jet.pt)
    #print "MET = %s with this phi = %s" %(met,met_phi)
    #for j, bjet in enumerate(bjets,1):
    #  print "bjet.pT = %s" %(bjets[0].pt)
    
    isGenuineTau = getIsGenuineTau(self,tau)
    #print "Tau has pt=%s, eta=%s, dz = %s, genPartFlav = %s"%(tau.pt, tau.eta, tau.dz, tau.genPartFlav)
    
    # Event Vars
    HT = 0.0
    LT = 0.0
    ST = 0.0
  
    for jet in jets: HT += jet.pt
    for e in electrons: LT += e.pt
    LT += tau.pt
    ST = HT + LT + MET
    
    #DiLepton Vars
    isOnZ = False;
    dileptonSystem_p4 = electron0.p4()+electron1.p4()
    dilepton_pt   = dileptonSystem_p4.Pt()
    dilepton_mass = dileptonSystem_p4.M()
    if( 75.0 < dilepton_mass < 105.0) : isOnZ = True
    
    ## check, jet to tau fake rate method
    IsOnZ             = isOnZ
    TauIsGenuine      = isGenuineTau
    TauPt             = tau.pt
    TauEta            = tau.eta
    TauDM             = tau.decayMode
    JetN              = len(jets)
    BJetN             = len(bjets)
    # HT                
    LeptonOnePt       = electron0.pt
    LeptonTwoPt       = electron1.pt
    DileptonPt        = dilepton_pt
    DileptonMass      = dilepton_mass
    DileptonDeltaEta  = dileptonSystem_p4.Eta()
    DileptonDeltaPhi  = deltaPhi(electron0.phi,electron1.phi)
    DileptonDeltaR    = electron0.DeltaR(electron1)
    
    self.IsOnZ[0]              = IsOnZ
    self.TauIsGenuine[0]       = TauIsGenuine
    self.TauPt[0]              = TauPt           
    self.TauEta[0]             = TauEta          
    self.TauDM[0]              = TauDM
    self.JetN[0]               = JetN            
    self.BJetN[0]              = BJetN           
    self.HT[0]                 = HT            
    self.LeptonOnePt[0]        = LeptonOnePt     
    self.LeptonTwoPt[0]        = LeptonTwoPt     
    self.DileptonPt[0]         = DileptonPt      
    self.DileptonMass[0]       = DileptonMass    
    self.DileptonDeltaEta[0]   = DileptonDeltaEta
    self.DileptonDeltaPhi[0]   = DileptonDeltaPhi
    self.DileptonDeltaR[0]     = DileptonDeltaR  
    
    #Fill Histos
    self.tree.Fill()
    
    return True

## Functions
def fillJetBranches(self,event,tau):
    """Help function to select jets and b tags, after removing overlap with tau decay candidates,
    and fill the jet variable branches."""
    
    #met     = self.met(event)
    jets,   bjets  = [ ], [ ]
    
    # SELECT JET, remove overlap with selected tau
    for jet in Collection(event,'Jet'):
      if abs(jet.eta)>2.4: continue #4.7: continue
      if jet.DeltaR(tau)<0.5: continue
      if jet.jetId<2: continue # Tight
      if jet.pt<self.jetCutPt: continue
      jets.append(jet)

      # B TAGGING
      if jet.btagDeepFlavB>self.deepjet_wp.medium and abs(jet.eta)<self.bjetCutEta:
        bjets.append(jet)
    
    # FILL JET BRANCHES
    #jets.sort( key=lambda j: self.ptnom(j),reverse=True)
    #bjets.sort(key=lambda j: self.ptnom(j),reverse=True)
    #self.out.njets[0]         = len(jets)
    #self.out.njets50[0]       = len([j for j in jets if self.ptnom(j)>50])
    #self.out.nfjets[0]        = nfjets
    #self.out.ncjets[0]        = ncjets
    #self.out.ncjets50[0]      = ncjets50
    #self.out.nbtag[0]         = nbtag
    
    # LEADING JET
    #if len(jets)>0:
    #  self.out.jpt_1[0]       = self.ptnom(jets[0])
    #  self.out.jeta_1[0]      = jets[0].eta
    #  self.out.jphi_1[0]      = jets[0].phi
    #  self.out.jdeepjet_1[0]  = jets[0].btagDeepFlavB
    #else:
    #  self.out.jpt_1[0]       = -1.
    #  self.out.jeta_1[0]      = -9.
    #  self.out.jphi_1[0]      = -9.
    #  self.out.jdeepjet_1[0]  = -9.
    
    # LEADING B JET
    #if len(bjets)>0:
    #  self.out.bpt_1[0]       = self.ptnom(bjets[0])
    #  self.out.beta_1[0]      = bjets[0].eta
    #else:
    #  self.out.bpt_1[0]       = -1.
    #  self.out.beta_1[0]      = -9.
    
    #return jets, bjets, met
    return jets, bjets

def getIsGenuineTau(self,tau):
  if not self.ismc : return False

  #KC checks
  #print "Hi there, you are in ModuleMuMuTau/getIsGenuineTau, under #KC checks-----"  
  #print "tau has genPartFlav = %s" %(tau.genPartFlav)
  if tau.genPartFlav==5: # genuine tau
    return True
  else:
    return False

#! /usr/bin/env python
# Author: Konstantinos Christoforou (Feb 2022)
# Description: Simple plotting script for jet-to-tau FakeRate studies
# ./plot_forTauFR.py -c mumutau --era UL2018
#from config.samples import *
from config.samples_forTauFR import *
from TauFW.Plotter.plot.string import filtervars
from TauFW.Plotter.plot.utils import LOG as PLOG

import TauFW.Plotter.tools.fakeFactors as fakeFactors


def plot(sampleset,channel,parallel=True,tag="",extratext="",outdir="plots",era="",
         varfilter=None,selfilter=None,fraction=False,pdf=False):
  """Test plotting of SampleSet class for data/MC comparison."""
  LOG.header("plot")
  
  # SELECTIONS
  if 'mumutau' in channel:
    #baseline = 'pt_mu0 > 10.0'
    #baseline = 'id_mu0 && id_mu1 && iso_mu0 < 0.15 && iso_mu1< 0.15 && iso_tau<0.15 && metfilter'
    #baseline = 'id_mu0 && id_mu1 && iso_mu0 < 0.15 && iso_mu1< 0.15 && id_tau >= 16 && iso_tau<0.15 && metfilter && JetN >= 2 && IsOnZ && BJetN >= 1'
    #baseline = 'id_mu0 && id_mu1 && iso_mu0 < 0.15 && iso_mu1< 0.15 && id_tau >= 16 && iso_tau<0.15 && metfilter && JetN >= 2 && IsOnZ && BJetN>=1'
    baseline = 'id_mu0 && id_mu1 && iso_mu0 < 0.15 && iso_mu1< 0.15 && id_tau >= 16 && iso_tau<0.15 && metfilter && JetN >= 2 && IsOnZ && BJetN>=1 && TauIsGenuine'
    #baseline = 'id_mu0 && id_mu1 && iso_mu0 < 0.15 && iso_mu1< 0.15 && id_tau >= 16 && iso_tau<0.15 && metfilter && JetN >= 2 && !IsOnZ && BJetN==0'
    #baseline = 'id_mu0 && id_mu1 && iso_mu0 < 0.15 && iso_mu1< 0.15 && id_tau >= 16 && iso_tau<0.15 && metfilter && JetN >= 2 && !IsOnZ && BJetN'
  elif 'mumu' in channel:
    baseline = "q_1*q_2<0 && pt_2>15 && iso_1<0.15 && iso_2<0.15 && idMedium_1 && idMedium_2 && !extraelec_veto && !extramuon_veto && metfilter && m_ll>20"
  elif 'mutau' in channel:
    baseline = "q_1*q_2<0 && iso_1<0.15 && idDecayModeNewDMs_2 && idDeepTau2017v2p1VSjet_2>=16 && idDeepTau2017v2p1VSe_2>=2 && idDeepTau2017v2p1VSmu_2>=8 && !lepton_vetoes_notau && metfilter"
  elif 'etau' in channel:
    baseline = "q_1*q_2<0 && iso_1<0.10 && mvaFall17noIso_WP90_1 && idDecayModeNewDMs_2 && idDeepTau2017v2p1VSjet_2>=16 && idDeepTau2017v2p1VSe_2>=8 && idDeepTau2017v2p1VSmu_2>=1 && !lepton_vetoes_notau && metfilter"
  elif 'tautau' in channel:
    baseline = "q_1*q_2<0 && iso_1<0.15 && iso_2<0.15 && idMedium_1 && idMedium_2 && !extraelec_veto && !extramuon_veto && m_ll>20 && metfilter"
  else:
    raise IOError("No baseline selection for channel %r defined!"%(channel))
  zttregion = "%s && mt_1<60 && dzeta>-25 && abs(deta_ll)<1.5"%(baseline) # && nbtag==0
  selections = [
    #Sel('baseline, no DeepTauVSjet',baseline.replace(" && idDeepTau2017v2p1VSjet_2>=16",""),only=["DeepTau"]),
    Sel("baseline",baseline),
    #Sel("baseline, pt > 50 GeV",baseline+" && pt_1>50"),
    #Sel("mt<60 GeV, dzeta>-25 GeV, |deta|<1.5",zttregion,fname="zttregion"),
    #Sel("0b",baseline+" && nbtag==0",weight="btagweight"),
    #Sel(">=1b",baseline+" && nbtag>=1",weight="btagweight"),
  ]
  
  #### DIFFERENTIAL pt/DM bins for TauID measurement
  ###wps = [
  ###  ('Medium',"idDeepTau2017v2p1VSjet_2>=16"),
  ###  #('VVVLoose && !Medium',"idDeepTau2017v2p1VSjet_2>=1 && idDeepTau2017v2p1VSjet_2<16")
  ###  #('0b, Medium',"idDeepTau2017v2p1VSjet_2>=16 && nbtag==0"),
  ###  #('0b, VVVLoose && !Medium',"idDeepTau2017v2p1VSjet_2>=1 && idDeepTau2017v2p1VSjet_2<16 && nbtag==0")
  ###  #('>=1b, Medium',"idDeepTau2017v2p1VSjet_2>=16 && nbtag>=1"),
  ###  #('>=1b, VVVLoose && !Medium',"idDeepTau2017v2p1VSjet_2>=1 && idDeepTau2017v2p1VSjet_2<16 && nbtag>=1")
  ###]
  ###pts = [20,30,40,50,70]
  ###dms = [0,1,10,11]
  ###for wp, wpcut in wps:
  ###  wpname  = wp.replace(" && !","-not")
  ###  basecut = baseline.replace("idDeepTau2017v2p1VSjet_2>=16",wpcut) #+" && nbtag==0"
  ###  for dm in dms:
  ###    name_ = "%s_dm%s"%(wpname,dm)
  ###    tit_  = "%s, DM%s"%(wp,dm)
  ###    cut_  = "%s && dm_2==%s"%(basecut,dm)
  ###    selections.append(Sel(name_,tit_,cut_,only=['m_vis','m_2'])) # DM bins
  ###  for i, ptlow in enumerate(pts):
  ###    if i<len(pts)-1: # ptlow < pt < ptup
  ###      ptup = pts[i+1]
  ###      name = "%s_pt%d-%d"%(wpname,ptlow,ptup)
  ###      tit  = "%s, %d < pt < %d GeV"%(wp,ptlow,ptup)
  ###      cut  = "%s && %s<pt_2 && pt_2<%s"%(basecut,ptlow,ptup)
  ###    else: # pt > ptlow (no upper pt cut)
  ###      name = "%s_pt%d-Inf"%(wpname,ptlow)
  ###      tit  = "%s, pt > %d GeV"%(wp,ptlow)
  ###      cut  = "%s && pt_2>%s"%(basecut,ptlow)
  ###    #selections.append(Sel(name,tit,cut,only=['m_vis','^m_2','mapRecoDM'])) # pt bins
  ###    for dm in dms:
  ###      name_ = "%s_dm%s"%(name,dm)
  ###      tit_  = "%s, DM%s"%(tit,dm)
  ###      cut_  = "%s && dm_2==%s"%(cut,dm)
  ###      selections.append(Sel(name_,tit_,cut_,only=['m_vis','^m_2'])) # pt-DM bins
  
  # VARIABLES
  variables = [
    #Var('TauPt'            , 'Tau pt'             , 60,  0   ,  300),
    Var('TauPt'            , 'Tau pt'             , [20, 25 , 30 , 40 , 50 , 80, 120]),
    Var('TauEta'           , 'Tau eta'            , 25, -2.5 ,   2.5),
    Var('TauDM'            , 'Tau DM'             , 12,  0   ,   12 ),
    Var('JetN'             , 'Number of jets'     , 10,  0   ,   10 ),
    Var('BJetN'            , 'Number of bjets'    ,  3,  0   ,    3 ),
    Var('HT'               , 'HT'                 , 30,  0   ,  600 ),
    Var('LT'               , 'LT'                 , 60,  0   ,  300 ),
    Var('ST'               , 'ST'                 , 50,  0   , 1000 ),
    Var('LeptonOnePt'      , 'Leading Muon pt'    , 50,  0   ,  150 , ctitle={'mumutau' : 'Leading Muon pt'} ),
    Var('LeptonTwoPt'      , 'subLeading Muon pt' , 50,  0   ,  100 ),
    Var('DileptonPt'       , 'Dilepton pt'        , 40,  0   ,  200 ),
    Var('DileptonMass'     , 'Dilepton mass'      , 20,  0   ,  200 ),
    Var('DileptonDeltaEta' , 'Dilepton deta'      , 23,  0   ,   4.6),
    Var('DileptonDeltaPhi' , 'Dilepton dphi'      , 32, -3.2 ,   3.2),
    Var('DileptonDeltaR'   , 'Dilepton dR'        , 23,  0   ,   4.6),

    #Var('pt_1',  "Muon pt",    40,  0, 120, ctitle={'etau':"Electron pt",'tautau':"Leading tau_h pt",'mumu':"Leading muon pt",'emu':"Electron pt"},cbins={"nbtag\w*>":(40,0,200)}),
    #Var('eta_1', "Muon eta",   30, -3,   3, ctitle={'etau':"Electron eta",'tautau':"Leading tau_h eta",'mumu':"Leading muon eta",'emu':"Electron eta"},ymargin=1.7,pos='T',ncols=2),
    #Var("jpt_1",  29,   10,  300, veto=[r"njets\w*==0"]),
    #Var("jeta_1", 53, -5.4,  5.2, ymargin=1.6,pos='T',ncols=2,veto=[r"njets\w*==0"]),
    #Var('njets',   8,  0,   8),
    #Var('nbtag', "Number of b jets (Medium, pt > 30 GeV)", 8, 0, 8),
    #Var('met',    50,  0, 150,cbins={"nbtag\w*>":(50,0,250)}),
    #Var('genmet', 50,  0, 150, fname="$VAR_log", logyrange=4, data=False, logy=True, ncols=2, pos='TT'),
    #Var('deta_ll', "deta(mutau_h)",  20, 0, 6.0, ctitle={'etau':"deta(etau_h)",'tautau':"deta(tautau)",'emu':"deta(emu)"},logy=True,pos='TRR',cbins={"abs(deta_ll)<":(10,0,3)}), #, ymargin=8, logyrange=2.6
  ]
  #if 'tau' in channel: # mutau, etau, tautau
  #  loadmacro("python/macros/mapDecayModes.C") # for mapRecoDM
  #  dmlabels  = ["h^{#pm}","h^{#pm}h^{0}","h^{#pm}h^{#mp}h^{#pm}","h^{#pm}h^{#mp}h^{#pm}h^{0}","Other"]
  #  variables += [
  #    Var('m_vis',          40,  0, 200, fname="mvis",ctitle={'mumu':"m_mumu",'emu':"m_emu"},cbins={"pt_\d>":(50,0,250),"nbtag\w*>":(60,0,300)},cpos={"pt_\d>[1678]0":'LL;y=0.88'}),
  #    Var("mapRecoDM(dm_2)", 5,  0,   5, fname="dm_2_label",title="Reconstructed tau_h decay mode",veto="dm_2==",position="TT",labels=dmlabels,ymargin=1.2),
  #  ]
  
  # PLOT
  selections = filtervars(selections,selfilter) # filter variable list with -V flag
  variables  = filtervars(variables,varfilter)  # filter variable list with -V flag
  outdir = ensuredir(repkey(outdir,CHANNEL=channel,ERA=era))
  exts   = ['png','pdf'] if pdf else ['png'] # extensions


  #for selection in selections:
  #  print ">>> Selection %r: %r"%(selection.title,selection.selection)
  #  #kc
  #  # stacks = sampleset.getstack(variables,selection,method='QCD_OSSS',parallel=parallel)
  #  stacks = sampleset.getstack(variables,selection,method='JetToTau_MisID',parallel=parallel)
  #  #stacks = sampleset.getstack(variables,selection,method='',parallel=parallel)
  #  fname  = "%s/$VAR_%s-%s-%s$TAG"%(outdir,channel.replace('mu','m').replace('tau','t'),selection.filename,era)
  #  text   = "%s: %s"%(channel.replace('mu',"#mu").replace('tau',"#tau_{h}"),selection.title)
  #  if extratext:
  #    text += ("" if '\n' in extratext[:3] else ", ") + extratext
  #  for stack, variable in stacks.iteritems():
  #    #position = "" #variable.position or 'topright'
  #    stack.draw(fraction=fraction)
  #    stack.drawlegend() #position)
  #    stack.drawtext(text)
  #    stack.saveas(fname,ext=exts,tag=tag)
  #    stack.close()

  selections_Fake = [
    Sel("baseline_Fake", baseline.replace("TauIsGenuine","!TauIsGenuine") ),
    ]
  selections_Fake = filtervars(selections_Fake,selfilter)
  selection_Fake  = selections_Fake[0]

  for selection in selections:
    print ">>> Selection %r: %r"%(selection.title,selection.selection)
    #kc
    #stacks = sampleset.getstack(variables,selection,method='QCD_OSSS',parallel=parallel)
    stacks = sampleset.getstack(variables,selection,method='JetToTau_MisID',parallel=parallel)    
    #stacks = sampleset.getstack(variables,selection,method='',parallel=parallel)    
    stacks_Fake = sampleset.getstack(variables,selection_Fake,method='',parallel=parallel)
    fname  = "%s/$VAR_%s-%s-%s$TAG"%(outdir,channel.replace('mu','m').replace('tau','t'),selection.filename,era)
    text   = "%s: %s"%(channel.replace('mu',"#mu").replace('tau',"#tau_{h}"),selection.title)
    if extratext:
      text += ("" if '\n' in extratext[:3] else ", ") + extratext
    for stack, variable in stacks.iteritems():
      stack.datahist.Reset()
      for stack_Fake, variable_Fake in stacks_Fake.iteritems():
        if variable_Fake == variable:
          stack.datahist.Add(stack_Fake.datahist)
      stack.draw(fraction=fraction)
      stack.drawlegend() #position)
      stack.drawtext(text)
      stack.saveas(fname,ext=exts,tag=tag)
      stack.close()



  '''
  # 17Feb22 checks, to be removed ///////////////////////////////////////////////////////////////////////////////////////////
   selections_0b_LnotTFake_forFRbins = []
  FRkeylist = []
  maxPt = 120
  etas = ["Barrel","Endcap"]
  prongs = ["1prong","3prong"]
  ptList = [20, 25 , 30 , 40 , 50 , 80]
   ## pt bins       
  dirList = ["plots/UL2018/TauFakeRate_UL2018_MuMuTau"]
  FRDict = {}
  ## prong bins
  for prong in prongs:
    name_ = "%s"%(prong)
    tit_  = "%s"%(prong)
    if prong == "1prong":
      cut_  = "%s && TauDM<5"%(baseline)
    else: ## prong == "3prong"                                                                       
      cut_  = "%s && TauDM>5"%(baseline)
    ## eta bins
    for eta in etas:
      name__ = "%s_eta%s"%(name_,eta)
      tit__  = "%s, %s"%(tit_,eta)
      if eta =="Barrel":
        cut__ = "%s && TauEta<1.5"%(cut_)
      else: # eta =="Endcap":                                                                 
        cut__ = "%s && TauEta<2.4"%(cut_)
      FakeRates = fakeFactors.FakeFactors(dirList, "Data", prong, eta, None, "mumutau")
      FRDictinPt  = {}
      FRDictinPt  = FakeRates.valuesDict
      firstvalue  = list(FRDictinPt.values())[0]
      for i, ptlow in enumerate(ptList):
        if i<len(ptList)-1: # ptlow < pt < ptup                                                       
          ptup = ptList[i+1]
          name___ = "%s_pt%d-%d"%(name__,ptlow,ptup)
          tit___  = "%s, %d < pt < %d GeV"%(tit__,ptlow,ptup)
          #cut_  = "%s && %s<TauPt && TauPt<%s"%(cuts_0b_LnotTFake,ptlow,ptup)
          cut___  = "%s && %s<TauPt && TauPt<%s"%(cut__,ptlow,ptup)
        else: # pt > ptlow (no upper pt cut)                                                        
          name___ = "%s_pt%d_%s"%(name__,ptlow,maxPt)
          tit___  = "%s, %d < pt < %d GeV"%(tit__,ptlow,maxPt)
          #cut_  = "%s && %s<TauPt && TauPt<%s"%(cuts_0b_LnotTFake,ptlow,maxPt)
          cut___  = "%s && %s<TauPt && TauPt<%s"%(cut__,ptlow,maxPt)
        selections_0b_LnotTFake_forFRbins.append(Sel(name___,tit___,cut___)) # pt-DM bins       
        FRDict[tit___] = firstvalue[i]
        FRkeylist.append(tit___)

  print FRDict
  #print list(FRDict.values())[0]
          
  for selection in selections_0b_LnotTFake_forFRbins:
    print selection.title
    FRDict_keyName = selection.title.replace("p_{T}","pt")
    print FRDict[FRDict_keyName]
    
  ## ================================================================================================= 
  # 17Feb22 checks, to be removed ///////////////////////////////////////////////////////////////////////////////////////////  
  '''
  

def main(args):
  channels  = args.channels
  eras      = args.eras
  parallel  = args.parallel
  varfilter = args.varfilter
  selfilter = args.selfilter
  notauidsf = args.notauidsf
  tag       = args.tag
  extratext = args.text
  fraction  = args.fraction
  pdf       = args.pdf
  outdir    = "plots/$ERA"
  fname     = "$PICODIR/$SAMPLE_$CHANNEL$TAG.root"
  for era in eras:
    for channel in channels:
      setera(era) # set era for plot style and lumi-xsec normalization
      addsfs = [ ] #"getTauIDSF(dm_2,genmatch_2)"]
      rmsfs  = [ ] if (channel=='mumu' or not notauidsf) else ['idweight_2','ltfweight_2'] # remove tau ID SFs
      split  = ['DY'] if 'tau' in channel else [ ] # split these backgrounds into tau components
      sampleset = getsampleset(channel,era,fname=fname,rmsf=rmsfs,addsf=addsfs,split=split) #,dyweight="")
      plot(sampleset,channel,parallel=parallel,tag=tag,extratext=extratext,outdir=outdir,era=era,
           varfilter=varfilter,selfilter=selfilter,fraction=fraction,pdf=pdf)
      sampleset.close()
  

if __name__ == "__main__":
  from argparse import ArgumentParser
  channels = ['mutau','etau','mumu','mumutau']
  eras = ['2016','2017','2018','UL2016_preVFP','UL2016_postVFP','UL2017','UL2018']
  description = """Simple plotting script for pico analysis tuples"""
  parser = ArgumentParser(prog="plot",description=description,epilog="Good luck!")
  parser.add_argument('-y', '--era',     dest='eras', nargs='*', choices=eras, default=['2017'],
                                         help="set era" )
  parser.add_argument('-c', '--channel', dest='channels', nargs='*', choices=channels, default=['mutau'],
                                         help="set channel" )
  parser.add_argument('-V', '--var',     dest='varfilter', nargs='+',
                                         help="only plot the variables passing this filter (glob patterns allowed)" )
  parser.add_argument('-S', '--sel',     dest='selfilter', nargs='+',
                                         help="only plot the selection passing this filter (glob patterns allowed)" )
  parser.add_argument('-s', '--serial',  dest='parallel', action='store_false',
                                         help="run Tree::MultiDraw serial instead of in parallel" )
  parser.add_argument('-F', '--fraction',dest='fraction', action='store_true',
                                         help="include fraction stack in ratio plot" )
  parser.add_argument('-p', '--pdf',     dest='pdf', action='store_true',
                                         help="create pdf version of each plot" )
  parser.add_argument('-r', '--nosf',    dest='notauidsf', action='store_true',
                                         help="remove DeepTau ID SF" )
  parser.add_argument('-t', '--tag',     default="", help="extra tag for output" )
  parser.add_argument('-T', '--text',    default="", help="extra text on plot" )
  parser.add_argument('-v', '--verbose', dest='verbosity', type=int, nargs='?', const=1, default=0, action='store',
                                         help="set verbosity" )
  args = parser.parse_args()
  LOG.verbosity = args.verbosity
  PLOG.verbosity = args.verbosity
  main(args)
  print "\n>>> Done."
  

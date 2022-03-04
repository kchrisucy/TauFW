#! /usr/bin/env python
# Author: Konstantinos Christoforou (Feb 2022)
# Description: Simple plotting script for jet-to-tau FakeRate studies
# ./plot_forTauFR.py -c mumutau --era UL2018
#from config.samples import *
from config.samples_forJetToTauFR import *
from TauFW.Plotter.plot.string import filtervars
from TauFW.Plotter.plot.utils import LOG as PLOG

import TauFW.Plotter.tools.fakeFactors as fakeFactors


def plot(sampleset,channel,parallel=True,tag="",extratext="",outdir="plots",era="",
         varfilter=None,selfilter=None,fraction=False,pdf=False,closure=False):
  """Test plotting of SampleSet class for data/MC comparison."""
  LOG.header("plot")
  
  # SELECTIONS
  if 'mumutau' in channel:
    if closure:
      #baseline = 'id_mu0 && id_mu1 && iso_mu0 < 0.15 && iso_mu1< 0.15 && (q_mu0*q_mu1) <= 0 && id_tau >= 16 && iso_tau<0.15 && metfilter && JetN >= 2 && IsOnZ && BJetN>=1 && TauIsGenuine'
      baseline = 'id_mu0 && id_mu1 && iso_mu0 < 0.15 && iso_mu1< 0.15 && (q_mu0*q_mu1) <= 0 && id_tau >= 16 && iso_tau<0.15 && metfilter && JetN >= 2 && IsOnZ && BJetN==0 && TauIsGenuine'
      #baseline = 'id_mu0 && id_mu1 && iso_mu0 < 0.15 && iso_mu1< 0.15 && (q_mu0*q_mu1) <= 0 && id_tau >= 16 && iso_tau<0.15 && metfilter && JetN >= 2 && IsOnZ && TauIsGenuine'
    else:
      #baseline = 'id_mu0 && id_mu1 && iso_mu0 < 0.15 && iso_mu1< 0.15 && (q_mu0*q_mu1) <= 0 && id_tau >= 1 && iso_tau<0.15 && metfilter && JetN >= 2 && IsOnZ && BJetN>=1'
      baseline = 'id_mu0 && id_mu1 && iso_mu0 < 0.15 && iso_mu1< 0.15 && (q_mu0*q_mu1) <= 0 && id_tau >= 1 && iso_tau<0.15 && metfilter && JetN >= 2 && IsOnZ && BJetN==0'
      #baseline = 'id_mu0 && id_mu1 && iso_mu0 < 0.15 && iso_mu1< 0.15 && (q_mu0*q_mu1) <= 0 && id_tau >= 1 && iso_tau<0.15 && metfilter && JetN >= 2 && IsOnZ'
  else:
    raise IOError("No baseline selection for channel %r defined!"%(channel))

  selections = [
    Sel("baseline",baseline),
  ]

  #### DIFFERENTIAL for TauID measurement #########
  if not closure :
    LooseTauWP = 'VVVLoose'
    TightTauWP = 'Medium'
    wps = [
      (LooseTauWP,'id_tau >= 1'),
      (TightTauWP,'id_tau >= 16'),
    ]

    for wp, wpcut in wps:
      basecut = baseline.replace("id_tau >= 1",wpcut)
      name = "%s"%(wp)
      tit  = "%s"%(wp)
      cut  = "%s"%(basecut)
      selections.append(Sel(name,tit,cut))
  
  # VARIABLES
  variables = [
    #Var('TauPt'            , 'Tau pt'             , 60,  0   ,  300),
    Var('TauPt'            , 'Tau pt'             , [20, 25 , 30 , 40 , 50 , 80, 120]),
    Var('TauEta'           , 'Tau eta'            , 25, -2.5 ,   2.5),
    Var('TauDM'            , 'Tau DM'             , 12,  0   ,   12 ),
    Var('JetN'             , 'Number of jets'     , 10,  0   ,   10 ),
    Var('BJetN'            , 'Number of bjets'    ,  3,  0   ,    3 ),
    Var('MET'              , 'MET'                , 60,  0   ,  300 ), ## in production
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
  ]
  
  # PLOT
  selections = filtervars(selections,selfilter) # filter variable list with -V flag
  variables  = filtervars(variables,varfilter)  # filter variable list with -V flag
  outdir = ensuredir(repkey(outdir,CHANNEL=channel,ERA=era))
  exts   = ['png','pdf'] if pdf else ['png'] # extensions

  ## for closure ccccccccccccccccccccccccc
  if closure:
    selections_Fake = [
      Sel("baseline_Fake", baseline.replace("TauIsGenuine","!TauIsGenuine") ),
    ]
    selections_Fake = filtervars(selections_Fake,selfilter)
    selection_Fake  = selections_Fake[0]
  ## ccccccccccccccccccccccccccccccccccccc

  for selection in selections:
    print ">>> Selection %r: %r"%(selection.title,selection.selection)
    ## ccccccccccccccccccccccccccccccccccccc
    if closure:
      stacks = sampleset.getstack(variables,selection,method='JetToTau_MisID',parallel=parallel)    
      stacks_Fake = sampleset.getstack(variables,selection_Fake,method='',parallel=parallel)
      ## ccccccccccccccccccccccccccccccccccc
    else:
      stacks = sampleset.getstack(variables,selection,method='',parallel=parallel)    
    fname  = "%s/$VAR_%s-%s-%s$TAG"%(outdir,channel.replace('mu','m').replace('tau','t'),selection.filename,era)
    text   = "%s: %s"%(channel.replace('mu',"#mu").replace('tau',"#tau_{h}"),selection.title)
    if extratext:
      text += ("" if '\n' in extratext[:3] else ", ") + extratext
    for stack, variable in stacks.iteritems():
      ## ccccccccccccccccccccccccccccccccccccc
      if closure:
        stack.datahist.Reset()
        for stack_Fake, variable_Fake in stacks_Fake.iteritems():
          if variable_Fake == variable:
            stack.datahist.Add(stack_Fake.datahist)
      ## ccccccccccccccccccccccccccccccccccccc
      stack.draw(fraction=fraction)
      stack.drawlegend() #position)
      stack.drawtext(text)
      stack.saveas(fname,ext=exts,tag=tag)
      stack.close()

## MAIN ####################################################
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
  closure   = args.closure
  if closure:
    outdir    = "plots/$ERA/Closure"
  else:
    outdir    = "plots/$ERA"
  fname     = "$PICODIR/$SAMPLE_$CHANNEL$TAG.root"
  for era in eras:
    for channel in channels:
      setera(era) # set era for plot style and lumi-xsec normalization
      addsfs = [ ] #"getTauIDSF(dm_2,genmatch_2)"]
      ## kc 23Feb22######
      rmsfs  = [ ]
      #rmsfs  = [ ] if (channel=='mumu' or not notauidsf) else ['idweight_2','ltfweight_2'] # remove tau ID SFs
      ###################
      split  = ['DY'] if 'tau' in channel else [ ] # split these backgrounds into tau components
      mergeMC = False
      sampleset = getsampleset(channel,era,fname=fname,rmsf=rmsfs,addsf=addsfs,split=split,mergeMC=mergeMC) #,dyweight="")
      plot(sampleset,channel,parallel=parallel,tag=tag,extratext=extratext,outdir=outdir,era=era,
           varfilter=varfilter,selfilter=selfilter,fraction=fraction,pdf=pdf,closure=closure)
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
  parser.add_argument('-C', '--closure', dest='closure', action='store_true', help="perform closure test" )
  parser.add_argument('-v', '--verbose', dest='verbosity', type=int, nargs='?', const=1, default=0, action='store',
                                         help="set verbosity" )
  args = parser.parse_args()
  LOG.verbosity = args.verbosity
  PLOG.verbosity = args.verbosity
  main(args)
  print "\n>>> Done."
  

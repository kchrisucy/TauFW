#! /usr/bin/env python
# Author: Konstantinos Christoforou (Feb 2022) 
# Description: Simple plotting script for pico analysis tuples
# ./plotJetToTauFRvariables_writeJson.py -c mumutau --era UL2018
#from config.samples import *
from config.samples_forJetToTauFR import *
from TauFW.Plotter.plot.string import filtervars
from TauFW.Plotter.plot.utils import LOG as PLOG
import TauFW.Plotter.tools.jsonWriter as jsonWriter 

#import os
import ROOT
import math
import ctypes 
import array

def plot(sampleset,channel,parallel=True,tag="",extratext="",outdir="plots",era="",
         varfilter=None,selfilter=None,fraction=False,pdf=False,verbose=0):
  """Test plotting of SampleSet class for data/MC comparison."""
  LOG.header("plot")
  
  # SELECTIONS
  if 'mumutau' in channel:
    baseline = 'id_mu0 && id_mu1 && iso_mu0 < 0.15 && iso_mu1< 0.15 && (q_mu0*q_mu1) <= 0 && id_tau >= 1 && iso_tau<0.15 && metfilter && JetN >= 2 && IsOnZ && BJetN==0'
    #baseline = 'id_mu0 && id_mu1 && iso_mu0 < 0.15 && iso_mu1< 0.15 && (q_mu0*q_mu1) <= 0 && id_tau >= 1 && iso_tau<0.15 && metfilter && JetN >= 2 && IsOnZ'
  else:
    raise IOError("No baseline selection for channel %r defined!"%(channel))

  selections = [
    Sel("baseline",baseline),
  ]
  
  #### DIFFERENTIAL pt/DM bins for TauID measurement
  LooseTauWP = 'VVVLoose'
  TightTauWP = 'Medium'
  LooseTauWPGenuine = LooseTauWP+'Genuine'
  TightTauWPGenuine = TightTauWP+'Genuine'
  wps = [
    (LooseTauWP,'id_tau >= 1'),
    (TightTauWP,'id_tau >= 16'),
    (LooseTauWPGenuine,'id_tau >= 1 && TauIsGenuine'),
    (TightTauWPGenuine,'id_tau >= 16 && TauIsGenuine'),
  ]

  #pts = [20,25,30,40,50,80,120]
  #pts = [20,30,80,120]
  maxPt = 120
  pts = [20]
  etas = ["Barrel","Endcap"]
  prongs = [1,3]
  ptList = [20, 25 , 30 , 40 , 50 , 80]
   
  for wp, wpcut in wps:
    basecut = baseline.replace("id_tau >= 1",wpcut ) #+" && nbtag==0"
    #print basecut
    for i, ptlow in enumerate(pts):
      if i<len(pts)-1: # ptlow < pt < ptup
        ptup = pts[i+1]
        name = "%s_pt%d-%d"%(wp,ptlow,ptup)
        tit  = "%s, %d < pt < %d GeV"%(wp,ptlow,ptup)
        cut  = "%s && %s<TauPt && TauPt<%s"%(basecut,ptlow,ptup)
      else: # pt > ptlow (no upper pt cut)
        name = "%s_pt%d-pt%s"%(wp,ptlow,maxPt)
        tit  = "%s, %d < pt < %d GeV"%(wp,ptlow,maxPt)
        cut  = "%s && %s<TauPt && TauPt<%s"%(basecut,ptlow,maxPt)
      #selections.append(Sel(name,tit,cut,only=['m_vis','^m_2','mapRecoDM'])) # pt bins
      for prong in prongs:
        name_ = "%s_prongs%s"%(name,prong)
        tit_  = "%s, prongs%s"%(tit,prong)
        if prong == 1:
          cut_  = "%s && TauDM<5"%(cut)
        else: ## prong == 3
          cut_  = "%s && TauDM>5"%(cut)
        #selections.append(Sel(name_,tit_,cut_,only=['TauPt'])) # pt-DM bins
        ###############################  
        for eta in etas:
          name__ = "%s_eta%s"%(name_,eta)
          tit__  = "%s, eta%s"%(tit_,eta)
          if eta =="Barrel":
            cut__ = "%s && abs(TauEta)<1.5"%(cut_)
          else: # eta =="Endcap":
            cut__ = "%s && abs(TauEta)>1.5 && abs(TauEta)<2.4"%(cut_)
          selections.append(Sel(name__,tit__,cut__,only=['TauPt'])) # pt-DM bins
  
  # VARIABLES
  variables = [
    Var('TauPt'            , 'Tau pt'             , 120, 20.0 ,  120.0),
    #Var('TauPt'            , 'Tau pt'             , [0, 20, 25 , 30 , 40 , 50 , 80, 120]),
    Var('TauEta'           , 'Tau eta'            , 25, -2.5 ,   2.5),
    Var('TauDM'            , 'Tau DM'             , 12,  0   ,   12 ),
  ]

  ## TauDM already 0,1,10 or 11 so if DM<5 -> 1 prong, if DM>5 -> 3 prong
  #hists_30_barrel_1prong = sample.gethist(variables,"TauPt<30 && TauEta<1.5 && TauDM<5") 
  
  # PLOT
  selections = filtervars(selections,selfilter) # filter variable list with -V flag
  variables  = filtervars(variables,varfilter)  # filter variable list with -V flag
  outdir = ensuredir(repkey(outdir,CHANNEL=channel,ERA=era))
  exts   = ['png','pdf'] if pdf else ['png'] # extensions
  LooseDataHists = []
  TightDataHists = []
  LooseMCHists = []
  TightMCHists = []
  LooseMCGenuineHists = []
  TightMCGenuineHists = []
  for selection in selections:
    #print ">>> Selection %r: %r"%(selection.title,selection.selection)
    stacks = sampleset.getstack(variables,selection,method='',parallel=parallel)
    fname  = "%s/$VAR_%s-%s-%s$TAG"%(outdir,channel.replace('mu','m').replace('tau','t'),selection.filename,era)
    text   = "%s: %s"%(channel.replace('mu',"#mu").replace('tau',"#tau_{h}"),selection.title)
    if extratext:
      text += ("" if '\n' in extratext[:3] else ", ") + extratext
    for stack, variable in stacks.iteritems():
      for h in stack.hists:
        if "prongs" not in selection.filename : continue
        if ("SingleMuon" not in h.GetName() and "Simulation" not in h.GetName()): continue 
        histo = h.Clone()
        #nBinsX = histo.GetNbinsX()
        histoName = histo.GetName()
        #print selection.filename
        #print histo.GetName()
        histoName += "_"+selection.filename
        histo.SetName(histoName)
        #histo = RebinHisto( histo  , [0.0]+[float(pt) for pt in ptList ]+[120.0])
        histo = RebinHisto( histo  , [float(pt) for pt in ptList ]+[120.0])
        #print histo.GetName()
        #print "h.GetName() =%s, nBinsX =%s "%(h.GetName(),nBinsX)
        if LooseTauWP in histoName:
          if "SingleMuon" in histoName : LooseDataHists.append(histo)
          elif "Simulation" in histoName : 
            if "Genuine" in histoName : LooseMCGenuineHists.append(histo)
            else :                      LooseMCHists.append(histo)
        elif TightTauWP in histoName:
          if "SingleMuon" in histoName : TightDataHists.append(histo)
          elif "Simulation" in histoName :
            if "Genuine" in histoName : TightMCGenuineHists.append(histo)
            else :                      TightMCHists.append(histo)
      ############
      stack.draw(fraction=fraction)
      stack.drawlegend() #position)
      stack.drawtext(text)
      stack.saveas(fname,ext=exts,tag=tag)
      stack.close()

  #for i in len(TightDataHists):
  #  numHistoData = TightDataHists[i]
  #  denHisto = LooseDataHists[i]    

  #print len(LooseDataHists), len(LooseMCHists), len(TightDataHists), len(TightMCHists)
  jsonList = []  
  jsonFileInfoDict = {}
  DM_eta_bins = []

  for prong in prongs:
    DM_eta_bin = "%sprong"%(prong)
    for eta in etas:
      #DM_eta_bin_ = "%s-%s"%(DM_eta_bin,eta)
      DM_eta_bin_ = "%s-%s"%(eta,DM_eta_bin)
      DM_eta_bins.append(DM_eta_bin_)

  for i, DM_eta_bin in enumerate(DM_eta_bins):
    numHistoData = TightDataHists[i]
    denHistoData = LooseDataHists[i]
    numHistoMC = TightMCHists[i]
    denHistoMC = LooseMCHists[i]
    numHistoGenuine = LooseMCGenuineHists[i]
    denHistoGenuine = TightMCGenuineHists[i]
    
    numHistoData.Add(numHistoGenuine,-1)
    denHistoData.Add(denHistoGenuine,-1)
    
    numHistoMC.Add(numHistoGenuine,-1)
    denHistoMC.Add(denHistoGenuine,-1)

    hEffData,jsonData  = GetEfficiencyHisto(numHistoData, denHistoData, "Data", outdir, era, channel, ptList, maxPt, verbose)
    hEffMC,jsonMC  = GetEfficiencyHisto(numHistoMC, denHistoMC, "Simulation", outdir, era, channel, ptList, maxPt, verbose)

    
    #myDelimiter = "-"
    dictKey     = channel + ":" + DM_eta_bin #+ myDelimiter.join(labelList[1:]) #exclude ptRange
    jsonFileInfoDict[dictKey + "-" + "Data"] = jsonData  
    jsonFileInfoDict[dictKey + "-" + "Simulation"] = jsonMC
    
  WriteJsonInfoToFile(jsonFileInfoDict,outdir,era,channel,verbose)
  



def RebinHisto(histo, newBinning):
  if isinstance(histo, ROOT.TH1):
    return histo.Rebin(len(newBinning)-1, histo.GetName(), array.array("d", newBinning))
    #return histo.Rebin(len(newBinning), histo.GetName(), array.array("d", newBinning))
  else:
    msg = "Must provide an object of type ROOT.TH1"
    raise Exception("!!!" + msg + "!!!")


def GetEfficiencyHisto(numHisto, denHisto, dsetType, outdir, era, channel, ptList, maxPt,verbose=0, statOpt=ROOT.TEfficiency.kFCP, confLevel=0.68): #confLevel=0.86638):
  '''
  the TEfficiency::SetStatisticOption() function sets the statistic option which affects the calculation of the confidence interval.

  statOpt (=statistics options)
  ==============================
  kFCP (=0)(default): using the Clopper-Pearson interval (recommended by PDG) sets kIsBayesian = false see also ClopperPearson
  kFNormal (=1)     : using the normal approximation sets kIsBayesian = false see also Normal
  kFWilson (=2)     : using the Wilson interval sets kIsBayesian = false see also Wilson
  kFAC (=3)         : using the Agresti-Coull interval sets kIsBayesian = false see also AgrestiCoull
  kFFC (=4)         : using the Feldman-Cousins frequentist method sets kIsBayesian = false see also FeldmanCousins
  kBJeffrey (=5)    : using the Jeffrey interval sets kIsBayesian = true, fBeta_alpha = 0.5 and fBeta_beta = 0.5 see also Bayesian
  kBUniform (=6)    : using a uniform prior sets kIsBayesian = true, fBeta_alpha = 1 and fBeta_beta = 1 see also Bayesian
  kBBayesian (=7)   : using a custom prior defined by fBeta_alpha and fBeta_beta sets kIsBayesian = true see also Bayesian
  kMidP (=8)        : using the Lancaster Mid-P method sets kIsBayesian = false

  confLevel (=confidence level)
  ==============================
  alpha = (1.0 - confLevel)/2
  2 sigma: confLevel 0.95 means alpha = 0.025
  1.64 sigma: confLevel 0.80 means alpha = 0.1
  1.5 sigma: confLevel 0.86638 means alpha = 0.06681
  1 sigma: confLevel 0.68 means alpha = 0.16
  
  '''
  # Definitions
  histoName = "Efficiency_%s_%s" % (dsetType, "b")
  #hEff = numHisto.Clone(histoName)
  hEff = numHisto.Clone()
  hEff.Reset()
  hEff.SetTitle(dsetType)
  hNum = ROOT.TH1D('Numerator'  ,'Numerator'  , 1, 0, 1)
  hDen = ROOT.TH1D('Denominator','Denominator', 1, 0, 1)
  
  # Construct info table (debugging)
  table  = []
  align  = "{:>6} {:>10} {:>10} {:>10} {:>12} {:>12} {:>12} {:^10} {:^10}"
  header = align.format("#", "pT-low", "pT-val", "pT-high", "Numerator", "Denominator", "Efficiency", "ErrorUp", "ErrorLow")
  hLine  = "="*100
  nBinsX = denHisto.GetNbinsX()
  table.append("{:^100}".format(histoName))
  table.append(hLine)
  table.append(header)
  table.append(hLine)

  # Define the dictionary for JSON file
  jsonInfoDelimiter = " & "
  jsonList    = []
  jsonList    = ["PtLowEdge__tight $\\tau_{h}$'s__loose $\\tau_{h}$'s__$\\varepsilon$__$\\sigma_{\\varepsilon}$__$\\sigma_{\\varepsilon}^{+}$__$\\sigma_{\\varepsilon}^{-}$__numerator__denominator"]
  jsonList[0] = jsonList[0].replace("__", jsonInfoDelimiter)


  # First determine the bin number [Bin# 0 contains the underflow. The Last bin (bin# nbins+1) contains the overflow.]
  #for j in range(1, nBinsX+1):
  for j in range(1, nBinsX+1):
    # Get the denominator
    denValue   = ctypes.c_double(0.0)
    denError_  = ctypes.c_double(0.0)
    denValue   = denHisto.IntegralAndError(j, j, denError_)
    denError   = denError_.value
    
    # Get the denominator
    numValue   = ctypes.c_double(0.0)
    numError_  = ctypes.c_double(0.0)
    numValue   = numHisto.IntegralAndError(j, j, numError_)
    numError   = numError_.value
    effValue   = 0.0
    effError   = 0.0
    histoName  = numHisto.GetName()
    histoTitle = numHisto.GetTitle()
    niceTitle  = histoTitle

    # Sanity
    if numValue < 0.0:
      numValueNew = 0.000001
      msg = "Integral of histogram %s for bin %d is less than zero (numValue = %.3f, denValue = %.3f). Setting numValue to %s" % (niceTitle, j, numValue, denValue, numValueNew)
      #self.Print(es + msg + ns, True)
      Print("!!!"+msg+"!!!", True)
      numValue = numValueNew
    if numError < 0.0:
      msg = "Error is less than zero (numError = %.3f, denError = %.3f)" % (numError, denError) 
      raise Exception("???" + msg + "???")
  
    # Numerator and Denominator histos
    if verbose:
      msg = "Setting bin content for numerator (denominator) to %.1f +/- %.1f (%.1f +/- %.1f)" % (numValue, numError, denValue, denError)
      Print("!!!"+msg+"!!!", True)
    
    # Sanity check (catch "nan")
    if (math.isnan(numValue) or math.isnan(numValue)):
      #aux.PrintTH1Info(numHisto)
      #aux.PrintTH1Info(denHisto)
      msg = "Setting bin content for numerator (denominator) to %.1f +/- %.1f (%.1f +/- %.1f)" % (numValue, numError, denValue, denError)
      raise Exception("???" + msg + "???")

    hNum.SetBinContent(1, numValue)
    hNum.SetBinError(1, numError)
    hDen.SetBinContent(1, denValue)
    hDen.SetBinError(1, denError)

    #print hNum.GetNbinsX()
    #print hDen.GetNbinsX()

    # Sanity check
    if not (ROOT.TEfficiency.CheckConsistency(hNum, hDen)):
      msg = "The input histograms are not valid / consistent for use with TEfficiency class"
      raise Exception("???" + msg + "???")

    # Calculate Efficiency (eff = Pass/Total) and error ( (eff*(1.0-eff))/sqrt(Denominator) )
    teff = ROOT.TEfficiency(hNum, hDen)

    # Set the Statistics options (https://root.cern.ch/doc/master/classTEfficiency.html)
    teff.SetStatisticOption(statOpt)
    teff.SetConfidenceLevel(confLevel) # 0.16=1 sigma corresponds to 68% CL, 0.025=2 sigma corresponds to 95% CL
    if verbose:
      msg = "teff.GetConfidenceLevel() = %s" %  (teff.GetConfidenceLevel() )
      Print("!!!"+msg+"!!!", True)
    
      # Get the efficiency & error (If the histograms are filled with weights, only bayesian methods and the normal approximation are supported)
    if verbose:
      msg = "TEfficiency confidence level fConfLevel=%.3f and the option fStatisticOption=%s"% (teff.GetConfidenceLevel(), teff.GetStatisticOption())
    effValue    = teff.GetEfficiency(1)
    effErrorUp  = teff.GetEfficiencyErrorUp(1)
    effErrorLow = teff.GetEfficiencyErrorLow(1)
    effError    = effErrorUp # legacy for backwards compatibility
    
    if numValue == 0 and denValue == 0:
      msg = "The efficiency calculation for bin %d is not possible (both numerator and denominator are zero). Efficiency = %.3f +/- %.3f." % (j, effValue, effError)
      Print("!!!"+msg+"!!!", True)
          
    # Bin-range or overflow bin?
    ptLowV = denHisto.GetXaxis().GetBinLowEdge(j)
    if ptLowV < 0:
      ptLowV = 0
    ptLow   = "%.1f" % (ptLowV)
    ptVal   = "%.1f" % (denHisto.GetXaxis().GetBinCenter(j))
    ptHigh  = "%.1f" % (denHisto.GetXaxis().GetBinUpEdge(j))
    ptRange = ptLow + " - " + ptHigh
    
    # Fill tefficiency object
    if math.isnan(effValue):
      msg = "Efficiency is nan! Forcing value to 0.0"
      Print("!!!"+msg+"!!!", True)
      effValue = 0.0
    hEff.SetBinContent(j, effValue)
    if math.isnan(effErrorUp):
      effErrorUp = 0.0
    if math.isnan(effErrorLow):
      effErrorLow = 0.0
    hEff.SetBinError(j, effErrorUp)

    # Save information in table
    row = align.format(j, ptLow, ptVal, ptHigh, "%.1f" % numValue, "%.1f" % denValue, "%.2f" % effValue, "%.2f" % effErrorUp, "%.2f" % effErrorLow, "%.2f" % numValue, "%.2f" % denValue)
    table.append(row)

    # For json file 
    json  = "%s__%.1f__%.1f__%.5f__%.5f__%.5f__%.5f__%.5f__%.5f" % (ptRange, numValue, denValue, effValue, effError, effErrorUp, effErrorLow, numValue, denValue)
    #print "json = ", json.replace("__", jsonInfoDelimiter)
    jsonList.append(json.replace("__", jsonInfoDelimiter))
      
    # Reset histos 
    hNum.Reset()
    hDen.Reset()
      
  # Finalise table
  table.append(hLine)

  # Print results
  #if verbose:
  if True:
    for i, line in enumerate(table):
      msg = line
      Print(msg, False)
      
  return hEff, jsonList


  # Write results to json file
  #  jsonFileInfoDict[dictKey + "-" + dsetType] = jsonList
  
def WriteJsonInfoToFile(jsonFileInfoDict,outdir,era,channel,verbose):
  outdirForJSON = "TauFakeRate_%s_%s" % (era, channel)
  outdir += "/" + outdirForJSON
  jsonWr = jsonWriter.JsonWriter(outdir, verbose)


  usedKeys = []
  for i, key in enumerate(jsonFileInfoDict.keys(), 1):
    infoList = jsonFileInfoDict[key]
    print infoList
    usedKey  = "-".join(key.split("-")[:-1])

    # For-loop: All info rows
    for k, row in enumerate(infoList, 1):
      if k == 1:
        continue
      else:
        param, value = _ConvertInfoToJsonParameter(key, row, k==2, k==(len(infoList)), usedKey in usedKeys)
        if verbose: Print("param = %s, value = %s" % (param, value), True)
        jsonWr.addParameter(param, value)
    usedKeys.append(usedKey)

  # Write the file
  #jsonFileName = "TauFakeRate_%s_%s" % (era, channel)
  jsonFileName = outdirForJSON + ".json"
  jsonWr.write(jsonFileName)
  jsonFilePath = outdir+jsonFileName
  return


def _ConvertInfoToJsonParameter(keyString, infoString, firstIndex=False, lastIndex=False, usedKey=False):

  jsonInfoDelimiter = " & "
  #self.Verbose("keyString = %s" % (keyString), True)
  print "keystring = ", keyString
  inputDir  = keyString.split(":")[0]
  print "inputDir = ", inputDir
  etaRegion = keyString.split(":")[-1].split("-")[0] #fixme: replace with a proper regex expression!
  print "etaRegion", etaRegion
  decayMode = keyString.split(":")[-1].split("-")[1] #decayMode = keyString.split("-")[1]
  print "decayMode" , decayMode
  
  print " keyString.split(:)[-1].split(-)" ,keyString.split(":")[-1].split("-")
  print " len keyString.split(:)[-1].split(-)" ,len( keyString.split(":")[-1].split("-") )

  if len(keyString.split(":")[-1].split("-")) > 3:
    #print " keyString.split(:)[-1].split(-)" ,keyString.split(":")[-1].split("-")
    ptRatio   = keyString.split("-")[2] 
    #print "ptRatio =  " , ptRatio
    dataType  = keyString.split(":")[-1].split("-")[3]
    #print "dataType = ", dataType
  else:
    ptRatio   = "N/A"
    #print "ptRatio =  " , ptRatio
    dataType  = keyString.split(":")[-1].split("-")[2]
    #print "dataType = ", dataType
    #self.Verbose("keyString.split(\"-\") = %s" %  (keyString.split("-")), True)
  
  #self.Verbose("InputDir = %s, DM = %s, etaRegion = %s, ptRatio = %s, dataType = %s" % (inputDir, decayMode, etaRegion, ptRatio, dataType), True)

  infoList  = infoString.split(jsonInfoDelimiter)
  ptRange     = infoList[0].replace(" - ", ", ")
  efficiency  = infoList[-6]
  error       = infoList[-5]
  errorUp     = infoList[-4]
  errorLow    = infoList[-3]
  numerator   = infoList[-2]
  denominator = infoList[-1]
  parameter   = "%s-%s_%s_%s" % (inputDir, decayMode, etaRegion, ptRatio)             
  parameter   = ConvertToFriendlySavename(parameter)
  parameter   = parameter.replace("_N/A", "") # in case ptRatio
  value       = ""
  #self.Verbose("PtRange = %s, Eff = %s, Error = %s, ErrorUp = %s, ErrorLow = %s" % (ptRange, efficiency, error, errorUp, errorLow), True)

  if firstIndex:
    if not usedKey: # needed because each pseudomulticrab has both "Data" and "Simulation"
      value += "{"
    value    += "\n\t\"%s\":{"         % (dataType)
  value    += "\n\t\t\"pt:[%s]\": {"     % (ptRange)
  value    += "\n\t\t\t\"value\": %s,"   % (efficiency)
  value    += "\n\t\t\t\"error\": %s,"   % (errorUp) # "error" can give nan -> crashes plotting script
  value    += "\n\t\t\t\"errorUp\": %s," % (errorUp)
  value    += "\n\t\t\t\"errorLow\": %s," % (errorLow)
  value    += "\n\t\t\t\"numerator\": %s," % (numerator)
  value    += "\n\t\t\t\"denominator\": %s" % (denominator)
  if lastIndex:
    value    += "\n\t\t}"
    value    += "\n\t}"
    if not usedKey:
      value+=","
  else:
    value += "\n\t\t},"
  return parameter, value


def ConvertToFriendlySavename(text):
  newText = text.replace("abs(TauJetEta)<1.5", "Barrel").replace("abs(TauJetEta)>1.5", "Endcap")
  newText = newText.replace("abs(TauJetEta)<1.46", "Barrel").replace("abs(TauJetEta)>1.46", "Endcap")
  newText = newText.replace("TauJetDM<5"   , "1prong")
  newText = newText.replace("TauJetDM=5..7", "2prong")
  newText = newText.replace("TauJetDM>7"   , "3prong")
  newText = newText.replace("TauJetPtRatio", "R").replace("=", "").replace("-", "_")
  newText = newText.replace("TauJetNProng<2", "1prong")
  newText = newText.replace("TauJetNProng>2", "3prong")
  return newText


## Nice Printing
#def _GetFName():         
#  fName = __file__.split("/")[-1].replace("pyc", "py")    
#  return fName
def Print(msg, printHeader=False):
  if printHeader==True:
    #print "=== ", _GetFName()
    print "=== ", __file__.split("/")[-1].replace("pyc", "py")    
    print "\t", msg
  else:
    print "\t", msg
  return
#################

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
  verbose   = args.verbosity
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
      mergeMC = True
      sampleset = getsampleset(channel,era,fname=fname,rmsf=rmsfs,addsf=addsfs,split=split,mergeMC=mergeMC) #,dyweight="")
      plot(sampleset,channel,parallel=parallel,tag=tag,extratext=extratext,outdir=outdir,era=era,
           varfilter=varfilter,selfilter=selfilter,fraction=fraction,pdf=pdf,verbose=verbose)
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
  

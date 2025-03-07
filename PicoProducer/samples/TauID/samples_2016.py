from TauFW.PicoProducer.storage.Sample import MC as M
from TauFW.PicoProducer.storage.Sample import Data as D
storage  = "/pnfs/psi.ch/cms/trivcat/store/user/ineuteli/samples/NANOAOD_2016/$PATH"
#storage  = "/pnfs/psi.ch/cms/trivcat/store/user/ineuteli/samples/nano/2016/$PATH"
store_T2 = "/pnfs/lcg.cscs.ch/cms/trivcat/store/user/areimers/NANOAOD/2016/LQTChannel/$PATH"
url      = None #"root://cms-xrd-global.cern.ch/"
filelist = "samples/files/2016/$SAMPLE.txt"
samples  = [
  
  # DRELL-YAN
  M('DY','DYJetsToLL_M-50',
    ###"/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16NanoAODv6-PUMoriond17_Nano25Oct2019_102X_mcRun2_asymptotic_v7_ext1-v1/NANOAODSIM", # OLD v6
    ###"/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16NanoAODv6-PUMoriond17_Nano25Oct2019_102X_mcRun2_asymptotic_v7_ext2-v1/NANOAODSIM", # OLD v6
    "/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16NanoAODv7-PUMoriond17_Nano02Apr2020_102X_mcRun2_asymptotic_v8_ext1-v1/NANOAODSIM",
    "/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16NanoAODv7-PUMoriond17_Nano02Apr2020_102X_mcRun2_asymptotic_v8_ext2-v1/NANOAODSIM",
    store=storage,url=url,files=filelist,opts='zpt=True',
  ),
  M('DY','DY1JetsToLL_M-50',
    ###"/DY1JetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16NanoAODv6-PUMoriond17_Nano25Oct2019_102X_mcRun2_asymptotic_v7-v1/NANOAODSIM", # OLD v6
    "/DY1JetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16NanoAODv7-PUMoriond17_Nano02Apr2020_102X_mcRun2_asymptotic_v8-v1/NANOAODSIM",
    store=storage,url=url,files=filelist,opts='zpt=True',
  ),
  M('DY','DY2JetsToLL_M-50',
    ###"/DY2JetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16NanoAODv6-PUMoriond17_Nano25Oct2019_102X_mcRun2_asymptotic_v7-v1/NANOAODSIM", # OLD v6
    "/DY3JetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16NanoAODv7-PUMoriond17_Nano02Apr2020_102X_mcRun2_asymptotic_v8-v1/NANOAODSIM",
    store=storage,url=url,files=filelist,opts='zpt=True',
  ),
  M('DY','DY3JetsToLL_M-50',
    ###"/DY3JetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16NanoAODv6-PUMoriond17_Nano25Oct2019_102X_mcRun2_asymptotic_v7-v1/NANOAODSIM", # OLD v6
    "/DY2JetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16NanoAODv7-PUMoriond17_Nano02Apr2020_102X_mcRun2_asymptotic_v8-v1/NANOAODSIM",
    store=storage,url=url,files=filelist,opts='zpt=True',
  ),
  M('DY','DY4JetsToLL_M-50',
    ###"/DY4JetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16NanoAODv6-PUMoriond17_Nano25Oct2019_102X_mcRun2_asymptotic_v7-v1/NANOAODSIM", # OLD v6
    "/DY4JetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16NanoAODv7-PUMoriond17_Nano02Apr2020_102X_mcRun2_asymptotic_v8-v1/NANOAODSIM",
    store=storage,url=url,files=filelist,opts='zpt=True',
  ),
  
  M('DY','DYJetsToLL_M-10to50',
    ###"/DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16NanoAODv6-PUMoriond17_Nano25Oct2019_102X_mcRun2_asymptotic_v7-v1/NANOAODSIM", # OLD v6
    "/DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16NanoAODv7-PUMoriond17_Nano02Apr2020_102X_mcRun2_asymptotic_v8-v1/NANOAODSIM",
    store=storage,url=url,files=filelist,opts=['zpt=True','useT1=False'],
  ),
  M('DY','DY1JetsToLL_M-10to50',
    ###"/DY1JetsToLL_M-10to50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16NanoAODv6-PUMoriond17_Nano25Oct2019_102X_mcRun2_asymptotic_v7-v1/NANOAODSIM", # OLD v6
    "/DY1JetsToLL_M-10to50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16NanoAODv7-PUMoriond17_Nano02Apr2020_102X_mcRun2_asymptotic_v8-v1/NANOAODSIM",
    store=storage,url=url,files=filelist,opts=['zpt=True','useT1=True'],
  ),
  M('DY','DY2JetsToLL_M-10to50',
    ###"/DY2JetsToLL_M-10to50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16NanoAODv6-PUMoriond17_Nano25Oct2019_102X_mcRun2_asymptotic_v7-v1/NANOAODSIM", # OLD v6
    "/DY2JetsToLL_M-10to50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16NanoAODv7-PUMoriond17_Nano02Apr2020_102X_mcRun2_asymptotic_v8-v1/NANOAODSIM",
    store=storage,url=url,files=filelist,opts=['zpt=True','useT1=True'],
  ),
  M('DY','DY3JetsToLL_M-10to50',
    ###"/DY3JetsToLL_M-10to50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16NanoAODv6-PUMoriond17_Nano25Oct2019_102X_mcRun2_asymptotic_v7-v1/NANOAODSIM", # OLD v6
    "/DY3JetsToLL_M-10to50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16NanoAODv7-PUMoriond17_Nano02Apr2020_102X_mcRun2_asymptotic_v8-v1/NANOAODSIM",
    store=storage,url=url,files=filelist,opts=['zpt=True','useT1=True'],
  ),
  M('DY','DY4JetsToLL_M-10to50',
    ###"/DY4JetsToLL_M-10to50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16NanoAODv6-PUMoriond17_Nano25Oct2019_102X_mcRun2_asymptotic_v7-v1/NANOAODSIM", # OLD v6
    "/DY4JetsToLL_M-10to50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16NanoAODv7-PUMoriond17_Nano02Apr2020_102X_mcRun2_asymptotic_v8-v1/NANOAODSIM",
    store=storage,url=url,files=filelist,opts=['zpt=True','useT1=True'],
  ),
  
  ## MASS-BINNED aMC@NLO
  #M('DY','DYJetsToLL_M-100to200',
  #  "/DYJetsToLL_M-100to200_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISummer16NanoAODv7-PUMoriond17_Nano02Apr2020_102X_mcRun2_asymptotic_v8_ext1-v1/NANOAODSIM",
  #  store=storage,url=url,files=None,opts=['zpt=True','useT1=True'],
  #),
  #M('DY','DYJetsToLL_M-200to400',
  #  "/DYJetsToLL_M-200to400_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISummer16NanoAODv7-PUMoriond17_Nano02Apr2020_102X_mcRun2_asymptotic_v8_ext2-v1/NANOAODSIM",
  #  store=storage,url=url,files=None,opts=['zpt=True','useT1=True'],
  #),
  #M('DY','DYJetsToLL_M-400to500',
  #  "/DYJetsToLL_M-400to500_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISummer16NanoAODv7-PUMoriond17_Nano02Apr2020_102X_mcRun2_asymptotic_v8_ext1-v1/NANOAODSIM",
  #  store=storage,url=url,files=None,opts=['zpt=True','useT1=True'],
  #),
  #M('DY','DYJetsToLL_M-500to700',
  #  "/DYJetsToLL_M-500to700_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISummer16NanoAODv7-PUMoriond17_Nano02Apr2020_102X_mcRun2_asymptotic_v8_ext1-v1/NANOAODSIM",
  #  store=storage,url=url,files=None,opts=['zpt=True','useT1=True'],
  #),
  #M('DY','DYJetsToLL_M-700to800',
  #  "/DYJetsToLL_M-700to800_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISummer16NanoAODv7-PUMoriond17_Nano02Apr2020_102X_mcRun2_asymptotic_v8_ext1-v1/NANOAODSIM",
  #  store=storage,url=url,files=None,opts=['zpt=True','useT1=True'],
  #),
  #M('DY','DYJetsToLL_M-800to1000',
  #  "/DYJetsToLL_M-800to1000_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISummer16NanoAODv7-PUMoriond17_Nano02Apr2020_102X_mcRun2_asymptotic_v8_ext1-v1/NANOAODSIM",
  #  store=storage,url=url,files=None,opts=['zpt=True','useT1=True'],
  #),
  #M('DY','DYJetsToLL_M-1000to1500',
  #  "/DYJetsToLL_M-1000to1500_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISummer16NanoAODv7-PUMoriond17_Nano02Apr2020_102X_mcRun2_asymptotic_v8_ext1-v1/NANOAODSIM",
  #  store=storage,url=url,files=None,opts=['zpt=True','useT1=True'],
  #),
  #M('DY','DYJetsToLL_M-1500to2000',
  #  "/DYJetsToLL_M-1500to2000_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISummer16NanoAODv7-PUMoriond17_Nano02Apr2020_102X_mcRun2_asymptotic_v8_ext1-v1/NANOAODSIM",
  #  store=storage,url=url,files=None,opts=['zpt=True','useT1=True'],
  #),
  #M('DY','DYJetsToLL_M-2000to3000',
  #  "/DYJetsToLL_M-2000to3000_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISummer16NanoAODv7-PUMoriond17_Nano02Apr2020_102X_mcRun2_asymptotic_v8_ext1-v1/NANOAODSIM",
  #  store=storage,url=url,files=None,opts=['zpt=True','useT1=True'],
  #),
  #M('DY','DYJetsToLL_M-3000toInf',
  #  "/DYJetsToLL_M-3000toInf_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIIFall17NanoAODv7-PU2017_12Apr2018_Nano02Apr2020_102X_mc2017_realistic_v8-v1/NANOAODSIM",
  #  store=None,url=None,files=None,opts=['zpt=True','useT1=True'],
  #),
  
  # HT-BINNED
  M('DY','DYJetsToLL_M-50_HT-70to100',
    "/DYJetsToLL_M-50_HT-70to100_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16NanoAODv7-PUMoriond17_Nano02Apr2020_102X_mcRun2_asymptotic_v8-v1/NANOAODSIM",
    store=storage,url=url,files=filelist,opts=['zpt=True','useT1=True'],
  ),
  M('DY','DYJetsToLL_M-50_HT-100to200',
    "/DYJetsToLL_M-50_HT-100to200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16NanoAODv7-PUMoriond17_Nano02Apr2020_102X_mcRun2_asymptotic_v8-v1/NANOAODSIM",
    "/DYJetsToLL_M-50_HT-100to200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16NanoAODv7-PUMoriond17_Nano02Apr2020_102X_mcRun2_asymptotic_v8_ext1-v1/NANOAODSIM",
    store=storage,url=url,files=filelist,opts=['zpt=True','useT1=True'],
  ),
  M('DY','DYJetsToLL_M-50_HT-200to400',
    "/DYJetsToLL_M-50_HT-200to400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16NanoAODv7-PUMoriond17_Nano02Apr2020_102X_mcRun2_asymptotic_v8-v1/NANOAODSIM",
    "/DYJetsToLL_M-50_HT-200to400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16NanoAODv7-PUMoriond17_Nano02Apr2020_102X_mcRun2_asymptotic_v8_ext1-v1/NANOAODSIM",
    store=storage,url=url,files=filelist,opts=['zpt=True','useT1=True'],
  ),
  M('DY','DYJetsToLL_M-50_HT-400to600',
    "/DYJetsToLL_M-50_HT-400to600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16NanoAODv7-PUMoriond17_Nano02Apr2020_102X_mcRun2_asymptotic_v8-v1/NANOAODSIM",
    "/DYJetsToLL_M-50_HT-400to600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16NanoAODv7-PUMoriond17_Nano02Apr2020_102X_mcRun2_asymptotic_v8_ext1-v1/NANOAODSIM",
    store=storage,url=url,files=filelist,opts=['zpt=True','useT1=True'],
  ),
  M('DY','DYJetsToLL_M-50_HT-600to800',
    "/DYJetsToLL_M-50_HT-600to800_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16NanoAODv7-PUMoriond17_Nano02Apr2020_102X_mcRun2_asymptotic_v8-v1/NANOAODSIM",
    store=storage,url=url,files=filelist,opts=['zpt=True','useT1=True'],
  ),
  M('DY','DYJetsToLL_M-50_HT-800to1200',
    "/DYJetsToLL_M-50_HT-800to1200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16NanoAODv7-PUMoriond17_Nano02Apr2020_102X_mcRun2_asymptotic_v8-v1/NANOAODSIM",
    store=storage,url=url,files=filelist,opts=['zpt=True','useT1=True'],
  ),
  M('DY','DYJetsToLL_M-50_HT-1200to2500',
    "/DYJetsToLL_M-50_HT-1200to2500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16NanoAODv7-PUMoriond17_Nano02Apr2020_102X_mcRun2_asymptotic_v8-v1/NANOAODSIM",
    store=storage,url=url,files=filelist,opts=['zpt=True','useT1=True'],
  ),
  M('DY','DYJetsToLL_M-50_HT-2500toInf',
    "/DYJetsToLL_M-50_HT-2500toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16NanoAODv7-PUMoriond17_Nano02Apr2020_102X_mcRun2_asymptotic_v8-v1/NANOAODSIM",
    store=storage,url=url,files=filelist,opts=['zpt=True','useT1=True'],
  ),
  
  # TTBAR
  M('TT','TT',
    ###"/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/RunIISummer16NanoAODv6-PUMoriond17_Nano25Oct2019_102X_mcRun2_asymptotic_v7-v2/NANOAODSIM", # OLD v6
    "/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/RunIISummer16NanoAODv7-PUMoriond17_Nano02Apr2020_102X_mcRun2_asymptotic_v8-v1/NANOAODSIM",
    store=storage,url=url,files=filelist,opts='toppt=True',
  ),
  
  # W+JETS
  M('WJ','WJetsToLNu',
    "/WJetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16NanoAODv6-PUMoriond17_Nano25Oct2019_102X_mcRun2_asymptotic_v7-v1/NANOAODSIM",
    "/WJetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16NanoAODv6-PUMoriond17_Nano25Oct2019_102X_mcRun2_asymptotic_v7_ext2-v1/NANOAODSIM",
    store=storage,url=url,files=filelist,
  ),
  M('WJ','W1JetsToLNu',
    "/W1JetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16NanoAODv6-PUMoriond17_Nano25Oct2019_102X_mcRun2_asymptotic_v7-v1/NANOAODSIM",
    store=storage,url=url,files=filelist,
  ),
  M('WJ','W2JetsToLNu',
    "/W2JetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16NanoAODv6-PUMoriond17_Nano25Oct2019_102X_mcRun2_asymptotic_v7-v1/NANOAODSIM",
    "/W2JetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16NanoAODv6-PUMoriond17_Nano25Oct2019_102X_mcRun2_asymptotic_v7_ext1-v1/NANOAODSIM",
    store=storage,url=url,files=filelist,
  ),
  M('WJ','W3JetsToLNu',
    "/W3JetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16NanoAODv6-PUMoriond17_Nano25Oct2019_102X_mcRun2_asymptotic_v7-v1/NANOAODSIM",
    "/W3JetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16NanoAODv6-PUMoriond17_Nano25Oct2019_102X_mcRun2_asymptotic_v7_ext1-v1/NANOAODSIM",
    store=storage,url=url,files=filelist,
  ),
  M('WJ','W4JetsToLNu',
    "/W4JetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16NanoAODv6-PUMoriond17_Nano25Oct2019_102X_mcRun2_asymptotic_v7-v1/NANOAODSIM",
    "/W4JetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16NanoAODv6-PUMoriond17_Nano25Oct2019_102X_mcRun2_asymptotic_v7_ext1-v1/NANOAODSIM",
    "/W4JetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16NanoAODv6-PUMoriond17_Nano25Oct2019_102X_mcRun2_asymptotic_v7_ext2-v1/NANOAODSIM",
    store=storage,url=url,files=filelist,
  ),
  
  # SINGLE TOP
  M('ST','ST_tW_antitop',
    "/ST_tW_antitop_5f_inclusiveDecays_TuneCP5_PSweights_13TeV-powheg-pythia8/RunIISummer16NanoAODv6-PUMoriond17_Nano25Oct2019_102X_mcRun2_asymptotic_v7-v1/NANOAODSIM",
    store=storage,url=url,files=filelist,
  ),
  M('ST','ST_tW_top',
    "/ST_tW_top_5f_inclusiveDecays_TuneCP5_PSweights_13TeV-powheg-pythia8/RunIISummer16NanoAODv6-PUMoriond17_Nano25Oct2019_102X_mcRun2_asymptotic_v7-v1/NANOAODSIM",
    store=storage,url=url,files=filelist,
  ),
  M('ST','ST_t-channel_antitop',
    "/ST_t-channel_antitop_4f_InclusiveDecays_TuneCP5_PSweights_13TeV-powheg-pythia8/RunIISummer16NanoAODv6-PUMoriond17_Nano25Oct2019_102X_mcRun2_asymptotic_v7-v1/NANOAODSIM",
    store=storage,url=url,files=filelist,
  ),
  M('ST','ST_t-channel_top',
    "/ST_t-channel_top_4f_InclusiveDecays_TuneCP5_PSweights_13TeV-powheg-pythia8/RunIISummer16NanoAODv6-PUMoriond17_Nano25Oct2019_102X_mcRun2_asymptotic_v7-v1/NANOAODSIM",
    store=storage,url=url,files=filelist,
  ),
  
  # DIBOSON
  M('VV','WW',
    "/WW_TuneCUETP8M1_13TeV-pythia8/RunIISummer16NanoAODv6-PUMoriond17_Nano25Oct2019_102X_mcRun2_asymptotic_v7-v1/NANOAODSIM",
    "/WW_TuneCUETP8M1_13TeV-pythia8/RunIISummer16NanoAODv6-PUMoriond17_Nano25Oct2019_102X_mcRun2_asymptotic_v7_ext1-v1/NANOAODSIM",
    store=storage,url=url,files=filelist,
  ),
  M('VV','WZ',
    "/WZ_TuneCUETP8M1_13TeV-pythia8/RunIISummer16NanoAODv6-PUMoriond17_Nano25Oct2019_102X_mcRun2_asymptotic_v7-v1/NANOAODSIM",
    "/WZ_TuneCUETP8M1_13TeV-pythia8/RunIISummer16NanoAODv6-PUMoriond17_Nano25Oct2019_102X_mcRun2_asymptotic_v7_ext1-v1/NANOAODSIM",
    store=storage,url=url,files=filelist,
  ),
  M('VV','ZZ',
    "/ZZ_TuneCUETP8M1_13TeV-pythia8/RunIISummer16NanoAODv6-PUMoriond17_Nano25Oct2019_102X_mcRun2_asymptotic_v7-v1/NANOAODSIM",
    "/ZZ_TuneCUETP8M1_13TeV-pythia8/RunIISummer16NanoAODv6-PUMoriond17_Nano25Oct2019_102X_mcRun2_asymptotic_v7_ext1-v1/NANOAODSIM",
    store=storage,url=url,files=filelist,
  ),
  
  # SINGLE MUON
  D('Data','SingleMuon_Run2016B',"/SingleMuon/Run2016B_ver2-Nano25Oct2019_ver2-v1/NANOAOD",
    store=storage,url=url,files=filelist,channels=["skim*",'mutau','mumu','emu'],
  ),
  D('Data','SingleMuon_Run2016C',"/SingleMuon/Run2016C-Nano25Oct2019-v1/NANOAOD",
    store=storage,url=url,files=filelist,channels=["skim*",'mutau','mumu','emu'],
  ),
  D('Data','SingleMuon_Run2016D',"/SingleMuon/Run2016D-Nano25Oct2019-v1/NANOAOD",
    store=storage,url=url,files=filelist,channels=["skim*",'mutau','mumu','emu'],
  ),
  D('Data','SingleMuon_Run2016E',"/SingleMuon/Run2016E-Nano25Oct2019-v1/NANOAOD",
    store=storage,url=url,files=filelist,channels=["skim*",'mutau','mumu','emu'],
  ),
  D('Data','SingleMuon_Run2016F',"/SingleMuon/Run2016F-Nano25Oct2019-v1/NANOAOD",
    store=storage,url=url,files=filelist,channels=["skim*",'mutau','mumu','emu'],
  ),
  D('Data','SingleMuon_Run2016G',"/SingleMuon/Run2016G-Nano25Oct2019-v1/NANOAOD",
    store=storage,url=url,files=filelist,channels=["skim*",'mutau','mumu','emu'],
  ),
  D('Data','SingleMuon_Run2016H',"/SingleMuon/Run2016H-Nano25Oct2019-v1/NANOAOD",
    store=storage,url=url,files=filelist,channels=["skim*",'mutau','mumu','emu'],
  ),
  
  # SINGLE ELECTRON
  D('Data','SingleElectron_Run2016B',"/SingleElectron/Run2016B_ver2-Nano25Oct2019_ver2-v1/NANOAOD",
    store=storage,url=url,files=filelist,channels=["skim*",'etau','ee'],
  ),
  D('Data','SingleElectron_Run2016C',"/SingleElectron/Run2016C-Nano25Oct2019-v1/NANOAOD",
    store=storage,url=url,files=filelist,channels=["skim*",'etau','ee'],
  ),
  D('Data','SingleElectron_Run2016D',"/SingleElectron/Run2016D-Nano25Oct2019-v1/NANOAOD",
    store=storage,url=url,files=filelist,channels=["skim*",'etau','ee'],
  ),
  D('Data','SingleElectron_Run2016E',"/SingleElectron/Run2016E-Nano25Oct2019-v1/NANOAOD",
    store=storage,url=url,files=filelist,channels=["skim*",'etau','ee'],
  ),
  D('Data','SingleElectron_Run2016F',"/SingleElectron/Run2016F-Nano25Oct2019-v1/NANOAOD",
    store=storage,url=url,files=filelist,channels=["skim*",'etau','ee'],
  ),
  D('Data','SingleElectron_Run2016G',"/SingleElectron/Run2016G-Nano25Oct2019-v1/NANOAOD",
    store=storage,url=url,files=filelist,channels=["skim*",'etau','ee'],
  ),
  D('Data','SingleElectron_Run2016H',"/SingleElectron/Run2016H-Nano25Oct2019-v1/NANOAOD",
    store=storage,url=url,files=filelist,channels=["skim*",'etau','ee'],
  ),
  
  # TAU
  D('Data','Tau_Run2016B',"/Tau/Run2016B_ver2-Nano25Oct2019_ver2-v1/NANOAOD",
    store=storage,url=url,files=filelist,channels=["skim*",'tautau'],
  ),
  D('Data','Tau_Run2016C',"/Tau/Run2016C-Nano25Oct2019-v1/NANOAOD",
    store=storage,url=url,files=filelist,channels=["skim*",'tautau'],
  ),
  D('Data','Tau_Run2016D',"/Tau/Run2016D-Nano25Oct2019-v1/NANOAOD",
    store=storage,url=url,files=filelist,channels=["skim*",'tautau'],
  ),
  D('Data','Tau_Run2016E',"/Tau/Run2016E-Nano25Oct2019-v1/NANOAOD",
    store=storage,url=url,files=filelist,channels=["skim*",'tautau'],
  ),
  D('Data','Tau_Run2016F',"/Tau/Run2016F-Nano25Oct2019-v1/NANOAOD",
    store=storage,url=url,files=filelist,channels=["skim*",'tautau'],
  ),
  D('Data','Tau_Run2016G',"/Tau/Run2016G-Nano25Oct2019-v1/NANOAOD",
    store=storage,url=url,files=filelist,channels=["skim*",'tautau'],
  ),
  D('Data','Tau_Run2016H',"/Tau/Run2016H-Nano25Oct2019-v1/NANOAOD",
    store=storage,url=url,files=filelist,channels=["skim*",'tautau'],
  ),
  
  ## LQ
  #M('LQ','SLQ_single_M1100_L1p0_old',
  #  "/LQ_Single_M1000/LegacyRun2_2016_deepTauIDv2p1/USER",
  #  store=storage,url=url,files=filelist,nfilesperjob=30,
  #),
  #M('LQ','SLQ_single_M1100_L1p0',
  #  "/SingleScalarLQ_InclusiveDecay_M-1100_L-1p0_TuneCUETP8M1_13TeV-madgraph-pythia8/RunIISummer16NanoAODv7-PUMoriond17_Nano02Apr2020_102X_mcRun2_asymptotic_v8-v1/NANOAODSIM",
  #  #store=storage,url=url,files=filelist,
  #),
  
  M('LQ','SLQ_nonres_M2500_Arne',
   "LQTChannelTauTau_HigherMasspoints_Scalar_M2500",
   store=store_T2,files=filelist,#url=url,
  ),
  M('LQ','SLQ_nonres_M3000_Arne',
   "LQTChannelTauTau_HigherMasspoints_Scalar_M3000",
   store=store_T2,files=filelist,#url=url,
  ),
  M('LQ','SLQ_nonres_M4000_Arne',
   "LQTChannelTauTau_HigherMasspoints_Scalar_M4000",
   store=store_T2,files=filelist,#url=url,
  ),
  M('LQ','SLQ_nonres_M5000_Arne',
   "LQTChannelTauTau_HigherMasspoints_Scalar_M5000",
   store=store_T2,files=filelist,#url=url,
   #blacklist=['NANOAOD_663.root','NANOAOD_624.root','NANOAOD_609.root']
  ),
  M('LQ','SLQ_nonres_M7000_Arne',
   "LQTChannelTauTau_HigherMasspoints_Scalar_M7000",
   store=store_T2,files=filelist,#url=url,
  ),
  M('LQ','SLQ_nonres_M10000_Arne',
   "LQTChannelTauTau_HigherMasspoints_Scalar_M10000",
   store=store_T2,files=filelist,#url=url,
  ),
  
  M('LQ','VLQ_nonres_M2500_Arne',
   "LQTChannelTauTau_HigherMasspoints_Vector_M2500",
   store=store_T2,files=filelist,#url=url,
  ),
  M('LQ','VLQ_nonres_M3000_Arne',
   "LQTChannelTauTau_HigherMasspoints_Vector_M3000",
   store=store_T2,files=filelist,#url=url,
  ),
  M('LQ','VLQ_nonres_M4000_Arne',
   "LQTChannelTauTau_HigherMasspoints_Vector_M4000",
   store=store_T2,files=filelist,#url=url,
  ),
  M('LQ','VLQ_nonres_M5000_Arne',
   "LQTChannelTauTau_HigherMasspoints_Vector_M5000",
   store=store_T2,files=filelist,#url=url,
  ),
  M('LQ','VLQ_nonres_M7000_Arne',
   "LQTChannelTauTau_HigherMasspoints_Vector_M7000",
   store=store_T2,files=filelist,#url=url,
  ),
  M('LQ','VLQ_nonres_M10000_Arne',
   "LQTChannelTauTau_HigherMasspoints_Vector_M10000",
   store=store_T2,files=filelist,#url=url,
  ),
  
]
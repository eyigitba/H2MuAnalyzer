import os
import sys
import array
from ROOT import *

class Plot_Config:

    def __init__(self, channel):
	self.channel = channel
    	self.signals = []
    	self.bkgs = []
    	self.data = ["Data"]
	self.colors = {}

	self.LoadColors()
	self.LoadSampleNames()

    def LoadSampleNames(self):
      self.signals.append("ttH")
      self.signals.append("ZH")
      self.signals.append("WH")
      self.signals.append("VBF")
      self.signals.append("ggH")

      if self.channel == "WH_3l":
	self.bkgs.append("others")
        self.bkgs.append("triboson")
        self.bkgs.append("tZq")
        self.bkgs.append("tW")
        self.bkgs.append("ttZ")
        self.bkgs.append("ttbar")
        self.bkgs.append("WW")
        self.bkgs.append("ZZ+ggZZ")
        self.bkgs.append("WZ")
        self.bkgs.append("DY")

      elif self.channel == "ZH_4l":
	self.bkgs.append("others")
        self.bkgs.append("triboson")
        self.bkgs.append("tZq")
        self.bkgs.append("ttZ")
        self.bkgs.append("ttbar")
        self.bkgs.append("WZ+WW")
        self.bkgs.append("ggZZ")
        self.bkgs.append("ZZ")
        self.bkgs.append("DY")

      else:
	print "channel is %s" %self.channel
	sys.exit()
    # end def LoadSampleNames(self):

    def LoadColors(self):
      if self.channel == "WH_3l":
	self.colors["Data"] = kBlack

        self.colors["ggH"] = kRed
        self.colors["VBF"] = kOrange + 7
        self.colors["ZH"] =  kBlue + 1
        self.colors["WH"] =  kGreen + 2
#        self.colors["WH_neg"] = kViolet + 1
        self.colors["ttH"]  = kPink + 6

	# trying to use blue for Z, green for W, yellow for top, red for g/q
        self.colors["DY"] =  kAzure + 7
        self.colors["ZZ+ggZZ"] =  kCyan - 7
        self.colors["WZ"] =  kGreen - 9
        self.colors["WW"] =  kSpring -1
        self.colors["ttbar"] = kYellow - 9
        self.colors["ttZ"] = kOrange - 9
        self.colors["tW"]  = kOrange + 6
        self.colors["tZq"] = kRed - 7
        self.colors["triboson"] = kViolet - 9
        self.colors["others"] = kPink + 6 

      elif self.channel == "ZH_4l":
	self.colors["Data"] = kBlack

        self.colors["ggH"] = kRed
        self.colors["VBF"] = kOrange + 7
        self.colors["ZH"] =  kBlue + 1
        self.colors["WH"] =  kGreen + 2
        self.colors["ttH"]  = kPink + 6

        self.colors["DY"] =  kAzure + 7
        self.colors["ZZ"] =  kCyan - 7
	self.colors["ggZZ"] =  kSpring -1
        self.colors["WZ+WW"] =  kGreen - 9
        self.colors["ttbar"] = kYellow - 9
        self.colors["ttZ"] = kOrange - 9
        self.colors["tZq"] = kRed - 7
        self.colors["triboson"] = kViolet - 9
        self.colors["others"] = kPink + 6 

      else:
        print "channel is %s" %self.channel
        sys.exit()
    # def LoadColors(self):




import numpy as np
import pandas as pd


class CellDischargeData:
    """
    Battery cell data from discharge test.
    """

    def __init__(self, path):
       
        data_frame = pd.read_csv(path, sep="\t")
        self.time = data_frame.iloc[:,0]
        self.voltage = data_frame.iloc[:,1]
        self.current = data_frame.iloc[:,2]

  
    @classmethod
    def process(cls, path):
        """
        Process the original discharge data for one section. This section of
        data is used for model development.
        """
        data = cls(path)

        data_frame = pd.read_csv(path, sep="\t")
          

        data.time  = data_frame.iloc[:,0]
        data.voltage = data_frame.iloc[:,1]
        data.current = data_frame.iloc[:,2]
      
        return data

    @classmethod
    def process_discharge_only(cls, path):
        """
        Process the original discharge data for just the discharge portion.
        """
        data = cls(path)
        data_frame = pd.read_csv(path, sep="\t")
          
        data.time  = data_frame.iloc[:,0]
        data.voltage = data_frame.iloc[:,1]
        data.current = data_frame.iloc[:,2]
                     
        return data

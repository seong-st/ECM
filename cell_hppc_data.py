import numpy as np
import pandas as pd


class CellHppcData:
    """
    Data from HPPC battey cell test.
    """

    def __init__(self, path, all_data=False):
        
        data_frame = pd.read_csv(path, sep="\t")

        self.time = data_frame.iloc[:,0]
        self.voltage = data_frame.iloc[:,1]
        self.current = data_frame.iloc[:,2]
        self.flags = data_frame.iloc[:,3]
               
    def get_indices(self):
        
        id_end_pulse = np.where(self.flags == 1)[0]
        id_start_rest = np.where(self.flags == 2)[0]
        id_end_rest = np.where(self.flags == 3)[0]
        id_start_pulse = np.where(self.flags == 4)[0]
        
        return id_end_pulse, id_start_rest, id_end_rest, id_start_pulse


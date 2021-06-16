# -*- coding: utf-8 -*-
"""
Created on Tue May 25 15:22:42 2021

@author: Alice Fermigier
"""

    def check_for_dropped_frames(self):
        """
            Checks if dropped frames were dropped during the recording
        """
        logging.info("Checking dropped frames...")
        camera_name = self.ca_video_tdms.split("(0)-")[-1].split(".tdms")[0]
        experiment_name = self.ca_video_tdms.split("(0)-")[0]
        notdropped = check_mantis_dropped_frames(self.data_folder, camera_name, experiment_name, 
                                skip_analog_inputs=True, verbose=False)
        if notdropped:
            logging.info("No frames were dropped")

        return camera_name, experiment_name
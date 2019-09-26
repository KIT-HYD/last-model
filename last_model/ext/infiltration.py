import numpy as np

from last_model._utils import ExtensionBase

from last_model.model_functions import (
    infilt_mtx,
    infilt_pfd
)

class Infiltration(ExtensionBase):
    identifier = 'Infiltration'
    
    def run(self):
        if self.last.qb_u > 0 or self.last.m_surface > 0:
            # Initialisation of input masses and surface storages
            self.last.m_surface = (np.abs(self.last.qb_u) * self.last.rho * 1 * self.last.dtc) + self.last.m_surface # rainwater initially stored in a virtual water surface storage
            self.last.m_slt_surface = self.last.m_slt_surface + (np.abs(self.last.qb_u) * 1 * self.last.dtc * self.last.Cw_eventparticles) # solute mass initially stored in a virtual solute surface storage
            self.last.Cw_eventparticles = self.last.m_slt_surface / (self.last.m_surface / self.last.rho)
 
            m_input1 = ((self.last.k[0] + self.last.ks) / 2) * ((np.abs(self.last.psi_init[0] - self.last.psi_init[1]) / self.last.dz) + 1) * 1 * self.last.rho * self.last.dtc # water mass input into soil matrix
            m_input2 = (self.last.pfd_k * (np.pi * self.last.pfd_r**2)) * self.last.rho * self.last.dtc * self.last.n_mak  # water mass input into pfd

            self.last.mass_totalwaterinput = self.last.mass_totalwaterinput + (np.abs(self.last.qb_u) * self.last.rho * 1 * self.last.dtc) # control parameter, sums up the total incoming water mass
            self.last.mass_totalsoluteinput = self.last.mass_totalsoluteinput + (np.abs(self.last.qb_u) * self.last.dtc * self.last.Cw_eventparticles) # control parameter, sums up the total incoming solute mass

            if self.last.Cw_eventparticles == 0: # concentration input after rainfall event dependent on masses in surface storages
                self.last.Cw_eventparticles = self.last.m_slt_surface / (self.last.m_surface / self.last.rho) 

            theta_count = 0 # scenario counter 

            # Scenario 1: Initially water and solutes flow into the soil matrix and residual overflow infiltrates simultaneously the pfd
            if np.round(self.last.theta[0],4) < (self.last.ths * 0.8):
     
                if m_input1 > self.last.m_surface:
                    m_input1 = self.last.m_surface
        
                if self.last.Cw_eventparticles * (m_input1 / self.last.rho) > self.last.m_slt_surface:
                    self.last.Cw_eventparticles = self.last.m_slt_surface / (m_input1 / self.last.rho)
          
                # calls function describing infiltration into soil matrix
                self.last.age_event, self.last.bla, self.last.Cw_event, self.last.m_slt_surface, self.last.m_surface, self.last.position_z_event, self.last.solute_test, self.last.theta, theta_count = infilt_mtx(t_mix_mean2=self.last.t_mix_mean2, t_mix_SD2=self.last.t_mix_SD2, t_mix_mean1=self.last.t_mix_mean1, t_mix_SD1=self.last.t_mix_SD1, t_mix_ratio=self.last.t_mix_ratio, age_event=self.last.age_event, bla=self.last.bla, Cw_event=self.last.Cw_event, Cw_eventparticles=self.last.Cw_eventparticles, dz=self.last.dz, m=self.last.m, m_input1=m_input1, m_slt_surface=self.last.m_slt_surface, m_surface=self.last.m_surface, position_z_event=self.last.position_z_event, rho=self.last.rho, solute_test=self.last.solute_test, theta=self.last.theta, time=self.last.time, z=self.last.z)
 
                # if the incoming precipitation mass is higher than the infiltrating mass(m_input1) it comes to a water accumulation of this residual water at the surface.This accumulation then flows immediately into the pfd with m_input2.
                if self.last.m_surface > 0:
     
                    if m_input2 > self.last.m_surface:
                        m_input2 = self.last.m_surface
     
                    # calls function describing infiltration into pfd
                    self.last.bla1, self.last.m_slt_surface, self.last.m_surface, self.last.pfd_position_z, self.last.pfd_theta, self.last.solute_test2, theta_count = infilt_pfd(bla1=self.last.bla1, Cw_eventparticles=self.last.Cw_eventparticles, m_input2=m_input2, m_slt_surface=self.last.m_slt_surface, m_surface=self.last.m_surface, n_mak=self.last.n_mak, pfd_dz=self.last.pfd_dz, pfd_m=self.last.pfd_m, pfd_position_z=self.last.pfd_position_z, pfd_r=self.last.pfd_r, pfd_theta=self.last.pfd_theta, pfd_z=self.last.pfd_z, rho=self.last.rho, solute_test2=self.last.solute_test2, theta_count=theta_count, time=self.last.time)

            # Scenario 2: If the top grid element of the soil matris is saturated the water and solutes infiltrate only the pfd
            if theta_count == 0 and np.round(self.last.theta[0],2) == self.last.ths and self.last.pfd_theta[0] < 1:
    
                if m_input2 > self.last.m_surface:
                    m_input2 = self.last.m_surface

                # calls function describing infiltration into pfd
                self.last.bla1, self.last.m_slt_surface, self.last.m_surface, self.last.pfd_position_z, self.last.pfd_theta, self.last.solute_test2, theta_count = infilt_pfd(bla1=self.last.bla1, Cw_eventparticles=self.last.Cw_eventparticles, m_input2=m_input2, m_slt_surface=self.last.m_slt_surface, m_surface=self.last.m_surface, n_mak=self.last.n_mak, pfd_dz=self.last.pfd_dz, pfd_m=self.last.pfd_m, pfd_position_z=self.last.pfd_position_z, pfd_r=self.last.pfd_r, pfd_theta=self.last.pfd_theta, pfd_z=self.last.pfd_z, rho=self.last.rho, solute_test2=self.last.solute_test2, theta_count=theta_count, time=self.last.time)


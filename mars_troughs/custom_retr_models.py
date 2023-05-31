#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Custom retreat models
# updated - do we even need splines?
"""
import numpy as np
from abc import abstractmethod
from mars_troughs.model import Model
from scipy.interpolate import InterpolatedUnivariateSpline as IUS

from mars_troughs.generic_model import (ConstantModel, 
                                        LinearModel,
                                        QuadModel, 
                                        CubicModel, 
                                        PowerLawModel)

class RetreatModel(Model):
    """
    Abstract class for retreat models, that have a method
    called :meth:`get_retreat_at_t` that returns retreat lag
    as a function of time.
    """

    prefix: str = "retr_"
    """All parameters of retreat models start with 'retreat'."""
    @abstractmethod
    def get_retreat_at_t(self, time: np.ndarray) -> np.ndarray:
        """
        Retreat as a function of time

        Args:
            time (np.ndarray): times at which we want to calculate the retreat.

        Output:
            np.ndarray of the same size as time input containing values of retreat.
        """
        #return self.eval(time)
        raise NotImplementedError # this coems from Acc model
   
class TimeDependentRetreatModel(RetreatModel):
    """
    A retreat rate model that depends on the time dependent parameter
    (likely solar insolation or obliquity), R(Var(t)).
    R is in m/year. Interpolated splines are created for the parameter as
    a function of time for faster integration.

    Args:
        times (np.ndarray): times at which the variable (solar insolation, obliquity) is known
                            (in years)
        parameter (np.ndarray): values of the time dependent variable
        (solar insolation (in W/m^2), obliquity (in degrees) )
    """

    def __init__(self, times: np.ndarray, variable: np.ndarray):
        self._times = times
        self._variable = variable
        self._var_data_spline = IUS(self._times, self._variable)
        self._int_var_data_spline = self._var_data_spline.antiderivative()
        self._var2_data_spline = IUS(self._times, self._variable ** 2)
        self._int_var2_data_spline = self._var2_data_spline.antiderivative()
        self._var3_data_spline = IUS(self._times, self._variable ** 3)
        self._int_var3_data_spline = self._var3_data_spline.antiderivative()

    def get_retreat_at_t(self, time: np.ndarray) -> np.ndarray:
        """
        Calculates the retreat rate at times "time".

        Args:
            time (np.ndarray): times at which we want to calculate R, in years.
        Output:
            np.ndarray of the same size as time input containing values of
            retreat rates R, in m/year

        """
        
        # what if: this is evaulated so it can be zero based on the negative times,
        # and then get_rt is set at
        
        re = self.eval(self._var_data_spline(time))
        #re = self.eval(self._var_data_spline(time))
        if np.any(re <0):
            mask = re < 0
            re[mask] = 0
        return re
        #    re_masked = np.zeros((np.size(re)))
        #    mask = re <0
        #    re_masked[mask] = 0
        #    re_masked[~mask] = re[~mask]
        #    return re_masked
        #else:
        #    return re

        #return self.eval(self._var_data_spline(time))
        

    # within acc model, retreat is set (get_xt), using retreat. 

class Constant_Retreat(TimeDependentRetreatModel, ConstantModel):
    """
    The retreat rate is constant and does not depend on time.

    Args:
        constant (float, optional): default is 1 millimeter. The lag
            thickness at all times.
    """

    def __init__(
        self,
        obl_times: np.ndarray,
        obliquity: np.ndarray,
        constant: float = 1e-6,
    ):
        super().__init__(obl_times, obliquity)  # note: `super` maps to the LagModel parent class
        ConstantModel.__init__(self, constant=constant)

    def get_rt(self, time: np.ndarray):
        """
        Calculates the retreat distance r (in m) traveled by a point
        in the center of the high side of the trough. This distance  is a
        function of the retreat rate A as r(t)=integral(AR(obl(t)), dt) or
        dy/dt=R(obl(t))

        Args:
            time (np.ndarray): times at which we want to calculate r, in years.
        Output:
            np.ndarray of the same size as time input containing values of
            the retreat distance r, in meters.

        """
        #re = (self.constant * time)

        #if any(re < 0):
        #    re_masked = np.zeros((np.size(re)))
        #    mask = re < 0
        #    re_masked[mask] = 0
        #    re_masked[~mask] = re[~mask]
        #    return re_masked
        #else:
        #    return re

        return (self.constant * time)
              

class Linear_Retreat(TimeDependentRetreatModel, LinearModel):
    def __init__(
        self,
        obl_times: np.ndarray,
        obliquity: np.ndarray,
        constant: float = 1e-6,#1e-6, # was e-6, then e-5
        slope: float = 1e-8, #1e-8, # was e-8, then 1e-6
    ):
        LinearModel.__init__(self, constant, slope)
        super().__init__(obl_times, obliquity)

    def get_rt(self, time: np.ndarray):
        """
        Calculates the retreat distance r (in m) traveled by a point
        in the center of the high side of the trough. This distance  is a
        function of the retreat rate R as r(t)=integral(R(obl(t)), dt) or
        dy/dt=R(obl(t))

        Args:
            time (np.ndarray): times at which we want to calculate r, in years.
        Output:
            np.ndarray of the same size as time input containing values of
            the retreat distance r, in meters.

        """
        
        re = self.eval(self._var_data_spline(time))

        if np.any(re <0):
            mask = re <= 0
            re[mask] = 0
            spline = IUS(self._times, re)
            int_var_spline = spline.antiderivative()
        else:
            int_var_spline = self._int_var_data_spline
            
        return (self.constant*time + (self.slope*(int_var_spline(time)-int_var_spline(0))))

       #     re_masked = np.zeros((np.size(re)))
       #     mask = re < 0
       #     re_masked[mask] = 0
       #     re_masked[~mask] = re[~mask]
       #     
       #     re_spline = IUS(self._times, re_masked)
       #     int_re_spline = re_spline.antiderivative()
       #     return (self.constant*time + (self.slope*(int_re_spline(time)-int_re_spline(0))))
        #else:
        #return (self.constant*time + (self.slope*(self._int_var_data_spline(time)-self._int_var_data_spline(0))))
            


        #return (self.constant*time + (self.slope*(int_re_spline(time)-int_re_spline(0))))
                                      #(re_int - self._int_var_data_spline(0))))
    
        #return -(
        #    self.constant * time
        #    + (
        #        self.slope
        #        * (self._int_var_data_spline(time) - self._int_var_data_spline(0))
        #    )
        #)
        
        #re = self._var_data_spline(time)

        #if np.any(self.eval(re) < 0):
        #    re_masked = np.zeros((np.size(re)))
        #    mask = self.eval(re) < 0
        #    re_masked[mask] = 0
        #    re_masked[~mask] = self.eval(re[~mask])
            #self._int_var_data_spline(time) = 0 #  = #re.antiderivative()
        
        #retreat_t = IUS(self._times, TimeDependentRetreatModel.get_retreat_at_t(self, time))
        #print(self._int_var_data_spline(0))
        #if np.any(((self._int_var_data_spline(time)- self._int_var_data_spline(0)))<0):
        #    print('negatives')
        #    spline_time = self._int_var_data_spline(time)
        #    locs = np.argwhere(spline_time  < 0)
        #    spline_time[locs] = 0
        
        
        # the issue is (time)- 0 gives delta Obl as negatives. 
        #re = (self.constant * time + (self.slope * (self._int_var_data_spline(time)- self._int_var_data_spline(0))))
                                      #(self._int_var_data_spline(time)- self._int_var_data_spline(0))))
                                      #(spline_time - self._int_var_data_spline(0) )))
              
           #(self._int_var_data_spline(time) - self._int_var_data_spline(0))) )
        #return re
    

class Quadratic_Retreat(TimeDependentRetreatModel,QuadModel):
    def __init__(
        self,
        obl_times: np.ndarray,
        obliquity: np.ndarray,
        constant: float = 1e-6,
        slope: float = 1e-8,
        quad: float = 1e-20, # was 1e-20
        ):
        
        QuadModel.__init__(self, constant, slope, quad)
        super().__init__(obl_times, obliquity)
        
    def get_rt(self, time: np.ndarray):
        """
        Calculates the vertical distance y (in m) at traveled by a point
        in the center of the high side of the trough. This distance  is a
        function of the accumulation rate A as y(t)=integral(A(ins(t)), dt) or
        dy/dt=A(ins(t))

        Args:
            time (np.ndarray): times at which we want to calculate y, in years.
        Output:
            np.ndarray of the same size as time input containing values of
            the vertical distance y, in meters.

        """
        output = (
            self.constant * time
            + (
                self.slope
                * (self._int_var_data_spline(time) - self._int_var_data_spline(0))
                
                + self.quad
                * (
                    self._int_var2_data_spline(time)
                    - self._int_var2_data_spline(0)
                )
            )
        )
        
        #mask = output < 0
        
        #output[mask] = 0
        
        return output
            
        """
        re = self.eval(self._var_data_spline(time))
        re2 = self.eval(self._var2_data_spline(time))

        if np.any(re <0):
            re_masked = np.zeros((np.size(re)))
            mask = re < 0
            re_masked[mask] = 0
            re_masked[~mask] = re[~mask]
            
            re_spline = IUS(self._times, re_masked)
            int_re_spline = re_spline.antiderivative()
            
            if np.any(re2 < 0):
                re2_masked = np.zeros((np.size(re2)))
                mask = re2 < 0
                re2_masked[mask] = 0
                re2_masked[~mask] = re2[~mask]
                
                re2_spline = IUS(self._times, re_masked**2)
                int_re2_spline = re2_spline.antiderivative()
                
                return (self.constant*time + (self.slope*(int_re_spline(time)-int_re_spline(0)))
                    + self.quad * (int_re2_spline(time) - int_re2_spline(0)
            ))

                
            else:
                return (self.constant*time + (self.slope*(int_re_spline(time)-int_re_spline(0)))
                    + self.quad * (self._int_var2_data_spline(time) - self._int_var2_data_spline(0)
            ))
        
        else:
            return (self.constant*time + 
                    (self.slope*(self._int_var_data_spline(time)-self._int_var_data_spline(0))
                     + self.quad* (self._int_var2_data_spline(time) - self._int_var2_data_spline(0)
            )))
          """  
    
    
    
    
        

class Cubic_Retreat(TimeDependentRetreatModel, CubicModel):
    def __init__(
        self,
        obl_times: np.ndarray,
        obliquity: np.ndarray,
        constant: float = 1e-6,
        slope: float = 1e-8,
        quad: float = 1e-20,
        cubic: float =1e-30,
        ):
        
        CubicModel.__init__(self, constant, slope, quad, cubic)
        super().__init__(obl_times, obliquity)
        
    def get_rt(self, time: np.ndarray):
        """
        Calculates the vertical distance y (in m) at traveled by a point
        in the center of the high side of the trough. This distance  is a
        function of the accumulation rate A as y(t)=integral(A(ins(t)), dt) or
        dy/dt=A(ins(t))

        Args:
            time (np.ndarray): times at which we want to calculate y, in years.
        Output:
            np.ndarray of the same size as time input containing values of
            the vertical distance y, in meters.

        """
        # ithink therse an issue here. 
        return (self.constant * time
                + 
                  (
                    self.slope
                    * (self._int_var_data_spline(time) 
                       - self._int_var_data_spline(0))
                    
                    + self.quad
                    * (self._int_var2_data_spline(time)
                        - self._int_var2_data_spline(0))
                    
                    + self.cubic
                    * (self._int_var3_data_spline(time)
                        - self._int_var3_data_spline(0))
                  )
               )

class PowerLaw_Retreat(TimeDependentRetreatModel, PowerLawModel):
    def __init__(
        self,
        obl_times: np.ndarray,
        obliquity: np.ndarray,
        coeff: float = 0.1,
        exponent: float = -2,
        ):
        
        PowerLawModel.__init__(self, coeff, exponent)
        super().__init__(obl_times, obliquity)
        

    def get_rt(self, time: np.ndarray):
        """
        Calculates the vertical distance y (in m) at traveled by a point
        in the center of the high side of the trough. This distance  is a
        function of the accumulation rate A as y(t)=integral(A(ins(t)), dt) or
        dy/dt=A(ins(t))

        Args:
            time (np.ndarray): times at which we want to calculate y, in years.
        Output:
            np.ndarray of the same size as time input containing values of
            the vertical distance y, in meters.

        """
        self._variable_exp = self._variable**self.exponent
        self._var_exp_data_spline = IUS(self._times, self._variable_exp )
        self._int_var_exp_data_spline = self._var_exp_data_spline.antiderivative()
        
        return (self.coeff*
                     (self._int_var_exp_data_spline(time)
                     -self._int_var_exp_data_spline(0)
                     )
                )
    
        
    
    
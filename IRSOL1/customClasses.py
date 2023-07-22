import enum

class eErrors (enum.Enum): 
    E_all_fine = 0
    E_arg_error = 1
    E_end_programm = 2
    E_unknown_error = 3


class roi:
    def __init__ (self,x,y,dx,dy):
        self.x = x
        self.y = y
        self.dx = dx
        self.dy = dy

class focus_param:
    def __init__ (self,speckles,sigma_noise,selected_frames,):
        speckles = speckles
        sigma_noise = sigma_noise
        selected_frames = selected_frames
 
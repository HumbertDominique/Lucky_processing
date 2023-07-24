import enum

class eErrors (enum.Enum): 
    E_all_fine = 1
    E_arg_error = 0
    E_end_programm = 2
    E_unknown_error = 3


class roi:
    def __init__ (self,x,y,dx,dy):
        self.x = x
        self.y = y
        self.dx = dx
        self.dy = dy

 
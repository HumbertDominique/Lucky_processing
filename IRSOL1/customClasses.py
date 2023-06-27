import enum

class programmState (enum.Enum):       
    startup = 0
    geometry_selection = 1
    geometry_cylinder_C = 2
    geometry_cylinder_A = 3
    save_file = 4
    close_file = 5
    Exit = 6
    
class eErrors (enum.Enum): 
    E_all_fine = 0
    E_arg_error = 1
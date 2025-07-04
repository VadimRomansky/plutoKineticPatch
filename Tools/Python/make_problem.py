import os
import sys
import shutil
import subprocess
import menu
import pluto_files_IO as pfIO
from define_problem import DefineProblem

class MakeProblem(object):
  def __init__(self, work_dir, pluto_dir, auto_update_def, auto_update_mkfl):
    """Create the makefile for the PLUTO code from the defintions header file.

    This class creates a makefile with all necessary information
    that would compile the code based on the problem defined by user
    in defintions.h file.
    In case the makefile is already present, this class will
    read those default values to re-create the file.

    **Inputs**:
      1. work_dir = Path to PLUTO code working directory
      2. pluto_dir = Path to PLUTO source directory
      3. auto-update_def  = Boolean that indicates auto-update of defintions.h.
      4. auto-update_mkfl = Boolean that indicates auto-update of makefile.

    **Output**:
      It generates a makefile for compilation.
    """
    self.work_dir = work_dir
    self.pluto_dir = pluto_dir
    self.auto_update = auto_update_mkfl
    self.mkfl_name = self.work_dir+'/makefile'        
    Dp = DefineProblem(self.work_dir, self.pluto_dir, auto_update_def)
    self.additional_files = Dp.additional_files
    self.additional_flags = Dp.additional_flags
    self.header_files = Dp.header_files
    self.pluto_path = Dp.pluto_path
    self.chomboflag = Dp.flag_dict['WITH-CHOMBO']
    self.high_order = Dp.flag_dict['WITH-HO']
    self.particles  = Dp.default[Dp.entries.index('PARTICLES')]

    try: 
      len(Dp.kromeoptstr)
    except AttributeError:
      pass
    else:
      krome_log = open('krome_config.out','w')
      os.chdir(Dp.krome_dir)
      kout = subprocess.Popen(Dp.kromeoptstr+' -unsafe', shell=True, stdout = krome_log)
      os.chdir(work_dir)
      self.row = 1
      menu.Print("> Krome Compilation Output written in [krome_config.out]", row=self.row,sleep=1)
      headgenstr = "python  "+self.pluto_dir+"/Src/Cooling/KROME/generate_cooling_header.py"
      hout = subprocess.Popen(headgenstr,shell=True, stdout=subprocess.PIPE)
      menu.Print("> Generated cooling.h in Src/Cooling/KROME", row=self.row+2,sleep=1)  

    self.SetArch()
    self.row = 1
    menu.Print ("> Generating makefile... ["+self.arch+"]", row=self.row,sleep=1)
    
    if self.chomboflag:
      self.ChomboMakeVars()
      self.makefile_template = '/Src/Templates/makefile.chombo'
    else:
      self.makefile_template = '/Src/Templates/makefile'

    self.UpdateMkflTemplate()
      
  
  def SetArch(self):
    """Sets the Architecture for compilation of code.

    This attribute of the MakeProblem class looks for the
    'ARCH' keyword in makefile and if not defined will define it
    based on user's choice of the Makefile configuration
    from the Config/ folder. If already defined then will use
    that architecture for compilation of Makefile.
    """
    mkfl_exits = os.path.exists(self.mkfl_name)
    if mkfl_exits:
      pf = pfIO.PlutoFiles(self.mkfl_name)
      scrh = pf.LocateString('ARCH')

    if self.auto_update == 0 or not mkfl_exits or len(scrh[0]) == 0:
      def_list = []
      entries  = os.listdir(self.pluto_dir + '/Config')
      for def_file in entries: 
        if (def_file.endswith('.defs')): def_list.append(def_file)

      def_list.sort()
      menu.SetTitle("Change makefile")
      self.arch = menu.Browse(def_list)  
      self.arch_string = 'ARCH         = '+ self.arch + '\n'
    else:                           
      self.arch_string = scrh[0][1]
      self.arch   = self.arch_string.split()[2]

  def ChomboMakeVars(self):
    """Adds CHOMBO specific vars in the Makefile.

    This method of the MakeProblem class does necessary
    modification to the makefile so as to accomodate
    compilation of chombo (AMR) specific variables.
    """
    pf = pfIO.PlutoFiles(self.work_dir+'/definitions.h')
    scrh = pf.LocateString('DIMENSIONS')
    dims = scrh[0][1].split()[2]
    chombo_config_string = 'DIM='+dims

    if '--with-chombo:' in sys.argv:
      i = sys.argv.index('--with-chombo:') + 1
      try:
        sys.argv[1:]
      except IndexError:
        print ("Additional Configration Details Required for '--with-chombo:' flag")
        sys.exit()
      else:
        for y in sys.argv[i:]:
          chombo_config_string += ' '+y

    self.row += 1
    menu.Print("  - Chombo config string: "+chombo_config_string,row=self.row) 
    self.row += 1
    menu.Print("  - creating make.vars...",row=self.row,sleep=0) 
    os.chdir(self.pluto_dir+"/Lib/Chombo-3.2/lib")
    os.system("make "+chombo_config_string+" vars > make.vars\n")
    os.system("cp make.vars "+self.work_dir+"\n")
    os.chdir(self.work_dir)
      
          
  def UpdateMkflTemplate(self):
    """
    Updates Makefile with additional flags, files and modular makefile paths.
    """
    shutil.copy(self.pluto_dir + self.makefile_template, self.mkfl_name)
    pf = pfIO.PlutoFiles(self.mkfl_name)
 
    # Set ARCH (architecture), PLUTO_DIR and VPATH
    pf.ReplaceWord('ARCH',self.arch_string, DelOld=True)
    pf.ReplaceWord('PLUTO_DIR','PLUTO_DIR    = '+self.pluto_dir+'\n',DelOld=True)

    if (not self.chomboflag):
      if self.high_order:
        scrh = pf.LocateString('SRC')
        ipos = scrh[0][0] + 1
        pf.InsertLine('SRC_HO       = $(SRC)/New/High_Order\n', ipos)
        vpath = ['./','$(SRC_HO)', '$(SRC)', '$(SRC)/Time_Stepping','$(SRC)/States']
      else:
        vpath = ['./','$(SRC)/New', '$(SRC)', '$(SRC)/Time_Stepping','$(SRC)/States']

      pf.ReplaceWord ("VPATH", "VPATH        = "+':'.join(vpath)+'\n', DelOld=True)

    # Insert additional CFLAGS
    scrh = pf.LocateString('Additional_CFLAGS_here')
    ipos = scrh[0][0] + 3
    for x in self.additional_flags:
      pf.InsertLine('CFLAGS += '+x+'\n', ipos)
      ipos = ipos + 1

    # Insert additional header and object files 
    scrh = pf.LocateString('Additional_header_files_here')
    ipos = scrh[0][0] + 3
    for x in self.header_files:
      pf.InsertLine('HEADERS += '+ x +'\n',ipos)
      ipos = ipos + 1

    scrh = pf.LocateString('Additional_object_files_here')
    ipos = scrh[0][0] + 3
    for x in self.additional_files:
      pf.InsertLine('OBJ += '+x + '\n', ipos)
      ipos = ipos + 1
    
    # Insert additional makefiles
    for x in self.pluto_path:
      pf.InsertLine('include $(SRC)/' + x + 'makefile' + '\n',ipos)
      ipos = ipos + 1

    # Add particle specific module: CR, LP or DUST
    if (self.particles == 'PARTICLES_LP'):
      pf.InsertLine('include $(SRC)/Particles/makefile_lp' + '\n',ipos)
      ipos = ipos + 1

    if (self.particles == 'PARTICLES_CR'):
      pf.InsertLine('include $(SRC)/Particles/makefile_cr' + '\n',ipos)
      ipos = ipos + 1

    if (self.particles == 'PARTICLES_DUST'):
      pf.InsertLine('include $(SRC)/Particles/makefile_dust' + '\n',ipos)
      ipos = ipos + 1
      
    if (self.particles == 'PARTICLES_MC'):
      pf.InsertLine('include $(SRC)/Particles/makefile_mc' + '\n',ipos)
      ipos = ipos + 1
    if (self.particles == 'PARTICLES_KIN'):
      pf.InsertLine('include $(SRC)/Particles/makefile_kin' + '\n',ipos)
      ipos = ipos + 1
    
#    for x in self.additional_flags:
#      pf.InsertLine('CFLAGS += '+x+'\n', ipos)
#      ipos = ipos + 1
       

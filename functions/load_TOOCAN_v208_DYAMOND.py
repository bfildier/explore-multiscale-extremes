### Importation des bibliotheques
import sys
import gzip

class MCS_IntParameters(object): 
  
    def __init__(self):
        self.DCS_number       		= 0
        self.INT_qltyDCS	 		= 0
        self.INT_classif			= 0
        self.INT_duration           = 0 
        self.INT_UTC_timeInit       = 0.
        self.INT_localtime_Init		= 0
        self.INT_lonInit			= 0
        self.INT_latInit			= 0
        self.INT_UTC_timeEnd		= 0
        self.INT_localtime_End		= 0
        self.INT_lonEnd				= 0
        self.INT_latEnd				= 0
        self.INT_velocityAvg		= 0		
        self.INT_distance			= 0
        self.INT_lonmin				= 0
        self.INT_latmin				= 0
        self.INT_lonmax				= 0
        self.INT_latmax				= 0
        self.INT_TbMin	            = 0
        self.INT_surfmaxPix_235K	= 0
        self.INT_surfmaxkm2_235K 	= 0 
        self.INT_surfmaxkm2_220K 	= 0
        self.INT_surfmaxkm2_210K  	= 0
        self.INT_surfmaxkm2_200K 	= 0
        self.INT_surfcumkm2_235K 	= 0 
        self.INT_classif_JIRAK 	    = 0
        self.INT_classif_MADDOX	    = 0

## Ajout B. Fildier ##

    def __repr__(self):
        """Creates a printable version of the Distribution object. Only prints the 
        attribute value when its string fits is small enough."""

        out = '< MCS_IntParameters object:\n'
        # print keys
        for k in self.__dict__.keys():
            out = out+' . %s: '%k
            if sys.getsizeof(getattr(self,k).__str__()) < 80:
                # show value
                out = out+'%s\n'%str(getattr(self,k))
            else:
                # show type
                out = out+'%s\n'%getattr(self,k).__class__
        out = out+' >'

        return out

## Fin ajout B. Fildier ##


class MCS_Lifecycle(object):

    def __init__(self):
        self.QCgeo_IRimage		= []
        self.LC_tbmin			= []
        self.LC_tbavg_235K 		= []
        self.LC_tbavg_208K      = []
        self.LC_tbavg_200K      = []
        self.LC_tb_90th         = []
        self.LC_UTC_time		= []
        self.LC_localtime		= []
        self.LC_lon				= []
        self.LC_lat				= []
        self.LC_x				= []
        self.LC_y				= []
        self.LC_velocity		= []
        self.LC_sminor_235K	    = []
        self.LC_smajor_235K	    = []
        self.LC_ecc_235K	    = []
        self.LC_orientation_235K= []
        self.LC_sminor_220K		= []
        self.LC_smajor_220K		= []
        self.LC_ecc_220K	    = []
        self.LC_orientation_220K	= []
        self.LC_surfPix_235K		= []
        self.LC_surfPix_210K		= []
        self.LC_surfkm2_235K		= []
        self.LC_surfkm2_220K		= []   
        self.LC_surfkm2_210K  		= []   
        self.LC_surfkm2_200K		= []   

## Ajout B. Fildier ##

    def __repr__(self):
        """Creates a printable version of the Distribution object. Only prints the 
        attribute value when its string fits is small enough."""

        out = '< MCS_Lifecycle object:\n'
        # print keys
        for k in self.__dict__.keys():
            out = out+' . %s: '%k
            if sys.getsizeof(getattr(self,k).__str__()) < 80:
                # show value
                out = out+'%s\n'%str(getattr(self,k))
            else:
                # show type
                out = out+'%s\n'%getattr(self,k).__class__
        out = out+' >'

        return out

## Fin ajout B. Fildier ##


def load_TOOCAN_DYAMOND(FileTOOCAN):

    lunit=gzip.open(FileTOOCAN,'rt')

    #
    # Read the Header
    ##########################
    header1 = lunit.readline()
    header2 = lunit.readline()
    header3 = lunit.readline()
    header4 = lunit.readline()
    header5 = lunit.readline()
    header6 = lunit.readline()
    header7 = lunit.readline()
    header8 = lunit.readline()
    header9 = lunit.readline()
    header10 = lunit.readline()
    header11 = lunit.readline()
    header12 = lunit.readline()
    header13 = lunit.readline()
    header14 = lunit.readline()
    header15 = lunit.readline()
    header16 = lunit.readline()
    header17 = lunit.readline()
    header18 = lunit.readline()
    header19 = lunit.readline()
    header20 = lunit.readline()
    header21 = lunit.readline()
    header22 = lunit.readline()
    header23 = lunit.readline()

    t=header11.split()
    temporalresolution = float(t[-2])


    data = []
    iMCS = -1
    lines = lunit.readlines()
    for iline in lines: 
        Values = iline.split()
        # print(iline,Values)
        # sys.exit()

        if(Values[0] == '==>'):

            #
			# Read the integrated parameters of the convective systems
			###########################################################
            data.append(MCS_IntParameters())
            iMCS = iMCS+1
            data[iMCS].DCS_number 			= int(Values[1])	    # Label of the convective system in the segmented images
            data[iMCS].INT_qltyDCS			= int(Values[2])	    # Quality control of the convective system 
            data[iMCS].INT_classif			= int(Values[3])	    # classif
            data[iMCS].INT_duration			= float(Values[4])    	# duration of the convective system (MCSMIP: duration in hours, )
            data[iMCS].INT_UTC_timeInit		= int(Values[5])		# time TU of initiation of the convective system
            data[iMCS].INT_localtime_Init	= int(Values[6])		# local time of inititiation
            data[iMCS].INT_lonInit			= float(Values[7])		# longitude of the center of mass at inititiation
            data[iMCS].INT_latInit			= float(Values[8])		# latitude of the center of mass at inititiation
            data[iMCS].INT_UTC_timeEnd		= int(Values[9])		# time TU of dissipation of the convective system
            data[iMCS].INT_localtime_End	= int(Values[10])		# local hour of dissipation
            data[iMCS].INT_lonEnd			= float(Values[11])		# longitude of the center of mass at dissipation
            data[iMCS].INT_latEnd			= float(Values[12])		# latitude of the center of mass at dissipation
            data[iMCS].INT_velocityAvg		= float(Values[13])		# average velocity during its life cycle(m/s)
            data[iMCS].INT_distance			= float(Values[14])		# distance covered by the convective system during its life cycle(km)
            data[iMCS].INT_lonmin			= float(Values[15])		# longitude min of the center of mass during its life cycle
            data[iMCS].INT_latmin           = float(Values[16])     # latitude min of the center of mass during its life cycle
            data[iMCS].INT_lonmax			= float(Values[17])		# longitude max of the center of mass during its life cycle
            data[iMCS].INT_latmax			= float(Values[18])		# latitude max of the center of mass during its life cycle
            data[iMCS].INT_TbMin			= float(Values[19])		# minimum Brigthness temperature (K)
            data[iMCS].INT_surfmaxPix_235K	= int(Values[20])		# maximum surface for a 235K threshold of the convective system during its life cycle (pixel)
            data[iMCS].INT_surfmaxkm2_235K	= float(Values[21])		# maximum surfacefor a 235K threshold of the convective system during its life cycle (km2)
            data[iMCS].INT_surfmaxkm2_220K	= float(Values[22])		# maximum surfacefor a 235K threshold of the convective system during its life cycle (km2)
            data[iMCS].INT_surfmaxkm2_210K	= float(Values[23])		# maximum surfacefor a 235K threshold of the convective system during its life cycle (km2)
            data[iMCS].INT_surfmaxkm2_200K	= float(Values[24])		# maximum surfacefor a 235K threshold of the convective system during its life cycle (km2)
            data[iMCS].INT_surfcumkm2_235K	= float(Values[25]) 	# integrated cumulated surface for a 235K threshold of the convective system during its life cycle (km2)		
            data[iMCS].INT_classif_JIRAK    = float(Values[26]) 	    # classif jirak
            data[iMCS].INT_classif_MADDOX   = float(Values[27]) 	    # classif jirak
            
            data[iMCS].clusters = MCS_Lifecycle()

            inc = 0
        else:
            #
            # Read the parameters of the convective systems 
            #along their life cycles
            ##################################################
            data[iMCS].clusters.QCgeo_IRimage.append(int(Values[0]))	    			# quality control on the Infrared image
            data[iMCS].clusters.LC_tbmin.append(float(Values[1]))	    			# min brightness temperature of the convective system at day TU (K)
            data[iMCS].clusters.LC_tbavg_235K.append(float(Values[2]))	    		# average brightness temperature of the convective system at day TU (K) 
            data[iMCS].clusters.LC_tbavg_208K.append(float(Values[3]))             # min brightness temperature of the convective system at day TU (K)
            data[iMCS].clusters.LC_tbavg_200K.append(float(Values[4]))	    	    # min brightness temperature of the convective system at day TU (K)
            data[iMCS].clusters.LC_tb_90th.append(float(Values[5]))	    		# min brightness temperature of the convective system at day TU (K)
            
            data[iMCS].clusters.LC_UTC_time.append(int(Values[6]))	    			# day TU 
            data[iMCS].clusters.LC_localtime.append(int(Values[7]))	    		# local hour (h)
            data[iMCS].clusters.LC_lon.append(float(Values[8]))	    			# longitude of the center of mass (°)
            data[iMCS].clusters.LC_lat.append(float(Values[9]))	    			# latitude of the center of mass (°)
            data[iMCS].clusters.LC_x.append(int(Values[10]))		    			# column of the center of mass (pixel)
            data[iMCS].clusters.LC_y.append(int(Values[11]))		    			# line of the center of mass(pixel)
            data[iMCS].clusters.LC_velocity.append(float(Values[12]))	    		# instantaneous velocity of the center of mass (m/s)
            
            data[iMCS].clusters.LC_sminor_235K.append(float(Values[13]))	    #
            data[iMCS].clusters.LC_smajor_235K.append(float(Values[14]))	    #
            data[iMCS].clusters.LC_ecc_235K.append(float(Values[15]))  	#
            data[iMCS].clusters.LC_orientation_235K.append(float(Values[16]))   	#
            data[iMCS].clusters.LC_sminor_220K.append(float(Values[17]))	    #
            data[iMCS].clusters.LC_smajor_220K.append(float(Values[18]))	    #
            data[iMCS].clusters.LC_ecc_220K.append(float(Values[19]))  	#
            data[iMCS].clusters.LC_orientation_220K.append(float(Values[20]))   	#
            data[iMCS].clusters.LC_surfPix_235K.append(int(Values[21]))	    	# surface of the convective system at time day TU (pixel)
            data[iMCS].clusters.LC_surfPix_210K.append(int(Values[22]))	    	# surface of the convective system at time day TU (pixel)
            data[iMCS].clusters.LC_surfkm2_235K.append(float(Values[23]))  		# surface of the convective system for a 235K threshold
            data[iMCS].clusters.LC_surfkm2_220K.append(float(Values[24]))  		# surface of the convective system for a 200K threshold
            data[iMCS].clusters.LC_surfkm2_210K.append(float(Values[25]))  		# surface of the convective system for a 210K threshold
            data[iMCS].clusters.LC_surfkm2_200K.append(float(Values[26]))  		# surface of the convective system for a 220K threshold

    lunit.close()
    return data    

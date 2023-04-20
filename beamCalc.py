import numpy as np
import matplotlib.pyplot as plt

def beamCalcI(beam):
    #Associating Local Variables with Beam Struct.
    webHeight = beam[0]
    webThickness = beam[1]
    flangeWidth = beam[2]
    flangeThickness = beam[3]
    webMat = beam[4]
    flangeMat = beam[5]

    #This code sets the parameters for glue shear. On dissimilar materials it takes the average. Otherwise it uses the glue shear for the relevant material.
    if(webMat != flangeMat):
        glueShearStrength = (1391+989)/2 #psi
    elif(webMat == 0): #Since the potential that the materials are dissimlar is already take care of we only need to check one.
        glueShearStrength = 1391 #psi
    else:
        glueShearStrength = 989 #psi

    #This code here sets the parameters for the web material.
    if(webMat == 0): #if the Web is made of Oak 
        webDensity = 0.0240 #lb/in^3
        webTensileStrength = 17873 #psi
        webModulus = 1800000 #psi
        webShearStrength = 2873 #psi
    else: #if the Web is made of Pine.
        webDensity = 0.0140 #psi
        webTensileStrength = 14327 #psi
        webModulus = 1500000 #psi
        webShearStrength = 1492 #psi

    #This code here sets the parameters for the flange material.
    if(flangeMat == 0): #if the Web is made of Oak 
        flangeDensity = 0.0240 #lb/in^3
        flangeTensileStrength = 17873 #psi
        flangeModulus = 1800000 #psi
        flangeShearStrength = 2873 #psi
    else: #if the Web is made of Pine.
        flangeDensity = 0.0140 #psi
        flangeTensileStrength = 14327 #psi
        flangeModulus = 1500000 #psi
        flangeShearStrength = 1492 #psi

    #Other parameters are set on this block.
    Length = 20 #in.

    ################# CALCULATING BEAM CONSTANTS ####################
    Mass = (Length*webHeight*webThickness)*webDensity + 2*(Length*flangeThickness*flangeWidth)*flangeDensity #Calculates the volume, multiplies by density for mass and sums each piece.

    if(webMat < flangeMat): #Handles dissimilar materials by converting each cross section into one of equivalent material.
        E_chosen = webModulus #If the web is made of oak and the flange is made of pine.
    elif(webMat > flangeMat):
        E_chosen = flangeModulus #If the web is made of pine and the flange is made of oak.
    else:
        E_chosen = flangeModulus
    n_web = E_chosen/webModulus #Finds scaling ratios.
    n_flange = E_chosen/flangeModulus #Scales the width.
    flangeWidth2 = flangeWidth/n_flange
    webThickness2 = webThickness/n_web
    I_zc =  (1/12)*((flangeWidth2*((flangeThickness*2+webHeight)**3)) - ((flangeWidth2-webThickness2)*(webHeight**3))) #in^4. Moment of inertia of the beam about the bending axis. Using the I_xc formula from 00B.
    Q_centroid = ((((webThickness))*(webHeight/2))*(webHeight/4))+((flangeThickness)*(flangeWidth)*(flangeThickness/2 + webHeight/2)) #Calculating first moment of area relative to the shear line.
    Q_joint = ((flangeThickness)*(flangeWidth)*(flangeThickness/2))
    ################# CALCULATING LOADING ####################
    #Flexural Strength and Tensile Strength can be presumed to be equal if the material is isotropic.
    #Wood is not, but since we were given tensile strength data measured along the grain of the wood and we will be making our beam along a straight grain we can presume they are the same.
    #The gist of what this code needs to do is find the failure loads for each possible loading case.
    #(Flexural, Wood Shear and Glue Shear) all have to be within 20% of the dominant failure mode

    Load = 0 #Start with the lower bound of the load.
    LoadMax = 5000000 #Upper bound of load testing.
    bendingFailureWeb = 0 #These are variables which store the failure points of each material if they occur.
    bendingFailureFlange = 0
    print("---Computing Bending Failure---")
    for x in range(Load,LoadMax):
        M_max = (0.4*x*12) #Again, fun quirk of the geometry. On the left hand side of the load, the shear force is 0.4x the load. Since M is integral of V, V_left * distance to load = M_max.
        BendingStress_Flange = (M_max/(I_zc*n_flange))*(((2*flangeThickness+webHeight))/2) #psi. Using the S = -(M/I*n_i)*c formula.
        BendingStress_Web = (M_max/(I_zc*n_web))*((webHeight)/2) #psi. Using the S = -(M/I*n_i)*c formula. Note that C has to be different since the maximum bending stress occurs at the junction between web and flange.
        if(BendingStress_Flange > flangeTensileStrength and bendingFailureFlange == 0): #Once a failure point is found (when the bending stress exceeds the tensile strength), store it in the variable. Continue the loop to search for other failure points.
            bendingFailureFlange = x
        if(BendingStress_Web > webTensileStrength and bendingFailureWeb == 0):
            bendingFailureWeb = x
    if(bendingFailureWeb < bendingFailureFlange): #Once out of the loop, we're going to look through each failure point to determine which occured first.
        print("Web fails first at " + str(bendingFailureWeb) + " lbf.")
        if(bendingFailureFlange != 0):
            print("Flange also fails at " + str(bendingFailureFlange) + " lbf.")
    elif(bendingFailureFlange < bendingFailureWeb):
        print("Flange fails first at " + str(bendingFailureFlange) + " lbf.")
        if(bendingFailureWeb != 0):
            print("Web also fails at " + str(bendingFailureWeb) + " lbf.")
    else:
        print("Bending failure does not occur within load range.")
    print("---Computing Shear Failure---")

    #At this point we're going repeat the same basic loop structure but calculating shear failure instead.
    shearFailureWeb = 0 #Same deal as above.
    shearFailureFlange = 0
    shearFailureGlue = 0
    for x in range(Load,LoadMax):
        V_max = -(0.6*x) #This is just because of the geometry of the situation. The reaction on the left hand side is 0.4x the load since the ratio of lengths between the supports is 1.5, and 0.4 = 0.6/1.5. The reaction of the right hand side (and thus V_max) is the remaining 0.6x the load.
        ShearStressWoodWeb = (-V_max*Q_centroid)/(n_web*I_zc*(webThickness/n_web)) #This is presuming the shear stress for the wood acts along the centroid.
        ShearStressWoodFlange = (-V_max*Q_joint)/(n_flange*I_zc*(flangeWidth/n_flange)) #This is presuming the shear stress for the flange pieces acts along the joint between them and the wood.
        if(ShearStressWoodWeb > webShearStrength and shearFailureWeb == 0): #Same deal as above, but we need to consider glue shear too.
            shearFailureWeb = x
        if(ShearStressWoodFlange > flangeShearStrength and shearFailureFlange == 0):
            shearFailureFlange = x
        if(ShearStressWoodFlange > glueShearStrength and shearFailureGlue == 0):
            shearFailureGlue = x
    if(shearFailureWeb != 0):
        print("Web fails in shear at " + str(shearFailureWeb) + " lbf.")
    else:
        print("Web does not fail in shear within test range.")
    if(shearFailureFlange != 0):
        print("Flange fails in shear at " + str(shearFailureFlange) + " lbf.")
    else:
        print("Flange does not fail in shear within test range.")
    if(shearFailureGlue !=0):
        print("Glue fails in shear at " + str(shearFailureGlue) + " lbf.")
    else:
        print("Glue does not fail in shear within test range.")

    ################# DOMINANT FAILURE MODE AND LOCATION ####################
    #Calculate where the beam fails, how it fails (which mode).
    print("---Computing Failure Mode and Delta---")
    stresses = {bendingFailureWeb: "Web Bending Failure", bendingFailureFlange: "Flange Bending Failure", shearFailureWeb: "Web Shear Failure", shearFailureFlange: "Flange Shear Failure", shearFailureGlue: "Glue Shear Failure"}
    stressesSorted = dict(sorted(stresses.items()))
    print("Dominant failure mode occurs first: " + str(stressesSorted))
    A = list(stressesSorted)[0]
    B = list(stressesSorted)[4]
    f_d = round(100*abs(A-B)/((A+B)/2),2)
    print("Largest Failure Delta is: " + str(f_d) +"%")
    if(f_d > 20):
        print("This is outside design range requirments, test new beam parameters please.")
    
    ################# Strength to Weight Ratio ####################

    print("---Computing Strength to Weight---")
    print("Beam Mass is: " + str(round(Mass,2)) + " lbs.")
    print("Strength to Weight Ratio is: " + str(round(A/Mass,2)) + ".")

    MaxDeflection = ((-A*8*12)/(6*E_chosen*I_zc*20))*(20^2-8^2-12^2)
    print("Max Deflection is: " + str(round(MaxDeflection,3)) + " in.")
    return 1

def paramCheck(beam):
    webHeight = beam[0]
    webThickness = beam[1]
    flangeWidth = beam[2]
    flangeThickness = beam[3]

    intCheck = 1 #This is the integrity check variable. If any issue is found the paramcheck function returns a 0 (check failed). Otherwise it returns a 1.

    if(webHeight or webThickness or flangeWidth or flangeThickness) < 0.1875: #This whole block of code checks each beam paramater to make sure its within the competition parameters.
        print("Some Material Dimension is under 3/16 in")
        intCheck = 0
    elif(webHeight + 2*flangeThickness) > 4:
        print("Beam Height exceeds 4 in maximum.")
        intCheck = 0
    elif(flangeWidth) > 2:
        print("Beam Width exceeds 2 in maximum.")
        intCheck = 0
    elif((webHeight + 2*flangeThickness)/flangeWidth) > 2:
        print("Beam exceeds maximum height to width ratio of 2.")
        intCheck = 0
    elif(webThickness or FlangeThickness) > 0.75:
        print("Beam cannot be made with 3/4 in stock material.")
        intCheck = 0
    elif((flangeWidth/flangeThickness) or (webHeight/webThickness)) > 8:
        print("Web or Flange Exceeds aspect Ratio limit of 8.")
        intCheck = 0
    else:
        print("Beam is within competition parameters!")
    return intCheck

def beamDesign(beam):
    print("---Beam Designer v0.1 by James Doss---")
    print("---Testing Beam Parameters---")
    passed = paramCheck(beam)
    if(passed):
        beamCalcI(beam)
    return 1

beamDesign([1.375,0.25,1.0625,0.3125,0,1])
#Feed in basic dimensions: [Web Height, Web Thickness, Flange Width, Flange Thickness, Web Material, Flange Material]
#All dims given in inches. Returns are in lbf. Material Options are 0 for Oak, 1 for Pine. Beam Options are 0 for a standard I beam and 1 for an inverted T beam.
#To run this program simply input your desired beam dimensions in the format given in the line directly above. If you are getting errors try checking the load simulation range (by default 0 lbs to 500000lbs).

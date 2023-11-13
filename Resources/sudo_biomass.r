#"sudo_biomass function is self defined function for biomass estiamtion based on the allometry functions from the "GlobAllomeTree" database#
#the website address is :http://www.globallometree.org #
#"InputDf" is one row of the allometry equation dataframe, class is dataframe#
#"SudoData" is a vector of pesudo DBH value, Unit is "cm"#

sudo_biomass = function(InputDf,SudoData,NameList)
{
    library(stringr)
    #every input data frame is one row of the allometry equation data frame#
    #every input data frame is one row of a data frame#
    if(is.data.frame(InputDf)&dim(InputDf)[1]==1)
    {
        print("Input data is appliable")
        if (InputDf$Unit_X=="cm")
        {
            SudoData            <- SudoData
            SudoData[SudoData<=InputDf$Min_X|SudoData>=InputDf$Max_X]  <- NA
        }else
        if (InputDf$Unit_X=="cm x 10")
        {
            SudoData            <- SudoData/10
            SudoData[SudoData<=InputDf$Min_X|SudoData>=InputDf$Max_X]  <- NA
        }else
        if (InputDf$Unit_X=="in")
        {
            SudoData            <- SudoData*0.3937
            SudoData[SudoData<=InputDf$Min_X|SudoData>=InputDf$Max_X]  <- NA
        }else
        if (InputDf$Unit_X=="in2")
        {
            #because the square is contatined in the equation, therefore, no transfroamtion to the square#
            SudoData            <- SudoData*0.3937
            SudoData[SudoData<=sqrt(InputDf$Min_X)|SudoData>=sqrt(InputDf$Max_X)]  <- NA
        }else
        if (InputDf$Unit_X=="m")
        {
            SudoData            <- SudoData/100
            SudoData[SudoData<=InputDf$Min_X|SudoData>=InputDf$Max_X]  <- NA
        }else
        if (InputDf$Unit_X=="mm")
        {
            SudoData            <- SudoData*100
            SudoData[SudoData<=InputDf$Min_X|SudoData>=InputDf$Max_X]  <- NA
        }else
        if (InputDf$Unit_X=="mm2")
        {
            #because the square is contatined in the equation, therefore, no transfroamtion to the square#
            SudoData            <- SudoData/100
            SudoData[SudoData<=sqrt(InputDf$Min_X)|SudoData>=sqrt(InputDf$Max_X)]  <- NA
        }
        dbh = SudoData
        dbh1 = SudoData
        c= SudoData
        print("Input data has been adapted to different Unit data")
        
        #access the equations in each row#
        per_equation             <- tolower(as.character(InputDf$Substitute_equation))
        #Find the "=" symbol in the equation string,the charater after "=" will be the start of the equation#
        start_position           <- regexpr("=",per_equation)[1]+1
        #Find the end point in the equation string#
        end_postion              <- nchar(per_equation)
        #Get the final equation for following calculations#
        calc_equ                 <- str_sub(per_equation, start = start_position, end = end_postion)
        
        #output type check and data transformation#
        #Transfer the NA output type to the orginal, and return the result directly#
        InputDf$Output_TR       <- as.character(InputDf$Output_TR)
        InputDf$Output_TR[is.na(InputDf$Output_TR)] <- "original"
        
        #print(calc_equ)
        
        if (InputDf$Output_TR=="log10")
        {
            print("log10")
            sud_biomass          <- 10^(eval(parse(text=calc_equ)))
        }else
        if (InputDf$Output_TR=="log"|InputDf$Output_TR=="Log")
        {
            print("ln")
            sud_biomass          <- exp(eval(parse(text=calc_equ)))
        }else
        if(InputDf$Output_TR=="log100")
        {
            print("log100")
            sud_biomass          <- 10^(eval(parse(text=calc_equ)))
        }else
        if(InputDf$Output_TR=="sqrt")
        {
            print("sqrt")
            sud_biomass          <- (eval(parse(text=calc_equ)))^2
        }else
        {
            print("original")
            sud_biomass          <- eval(parse(text=calc_equ))
            #sud_biomass         <- rep(0,500)
        }
        
        #For the output Units diverse from g to lb, we have to uniform them to kg#
        if (InputDf$Unit_Y=="g")
        {
            sud_biomass           <- sud_biomass/1000
        }else
        if(InputDf$Unit_Y=="kg")
        {
            sud_biomass           <- sud_biomass
        }else
        if(InputDf$Unit_Y=="lb")
        {
            sud_biomass           <- sud_biomass*0.453592
        }else
        if(InputDf$Unit_Y=="Mg")
        {
            sud_biomass           <- sud_biomass/1000000
        }
        #format the calculation to a named row of data frame#
        sud_biomass               <- t(as.data.frame(sud_biomass))
        sud_biomass_df            <- data.frame(ID=InputDf$ID_AE,LAT=InputDf$Latitude,LON=InputDf$Longitude,Ecoregion_WWF=InputDf$Ecoregion_WWF,InputDf$Substitute_equation,sud_biomass)
        names(sud_biomass_df)     <- c("ID","LAT","LON","Ecoregion_WWF","Equations",NameList)
        
        return(sud_biomass_df)
    }else
    {
        print("Input data is not appliable, please check the requirment of the input data type")
    }
}


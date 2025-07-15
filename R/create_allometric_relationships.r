create_allometric_relationships = function(p) {

    det = snowcrab.db( DS="det.initial"  )

    allometry.snowcrab( "cw.mass", "male", redo=T )
    allometry.snowcrab( "chela.mass", "male", redo=T  )
    allometry.snowcrab( "cw.chela.mat", "male", redo=T  )
    allometry.snowcrab( "cw.mass", "female", redo=T  )
    allometry.snowcrab( "chela.mass", "female", redo=T  )
    allometry.snowcrab( "cw.chela.mat", "female", redo=T  )


    #Identify morphology errors and print, save to CSV
    yr.errors <- p$year.assessment
    fn.errors = file.path(project.datadirectory("bio.snowcrab"), "output", "morphology.errorsrrors")
    dir.create(fn.errors, recursive=TRUE, showWarnings=F)
    outfile.errors =  file.path( fn.errors, paste("morphologyerrors", yr.errors, ".csv", sep=""))
    outfile.errors2 =  file.path( fn.errors, paste("morphologyerrors.allyears", yr.errors, ".csv", sep=""))

    #Sex.errors: Unknown Sex
    sex.errors <- det[which(det$sex==sex.unknown),]
    if ( !is.na(sex.errors$trip[1])) sex.errors$error <- 'sex.errors'

    #Cw.errors: Carapace Width below 5 or greater than 185
    cw.errors <- det[ which(det$cw<5 | det$cw>185 ),]
    if ( !is.na(cw.errors$trip[1])) cw.errors$error <- 'cw.errors'

    #Chela.errors: Chela less than 1 or greater than 50
    chela.errors <- det[which(det$chela < 1 | det$chela > 50  ),]
    if ( !is.na(chela.errors$trip[1])) chela.errors$error <- 'chela.errors'

    #Abdomen.errors:Abdomen less than 1 and greater than 66
    abdomen.errors <- det[which(det$abdomen < 1 | det$abdomen > 66 ),]

    #abdomen.errors$error <- 'abdomen.errors' #BZ 2018 no abdomen lengths met "error" condition, broke script #
    if ( !is.na(abdomen.errors$trip[1]))  abdomen.errors$error='abdomen.errors' #replaced above statement

    #Mass.errors: Mass less than 1 or greater than 1500
    mass.errors <- det[which( det$mass < 1 | det$mass > 1500  ),]
    if ( !is.na(mass.errors$trip[1]))  mass.errors$error <- 'mass.errors'

    #Sex.a: Indeterminate sex based on measurements taken (abdomen values where sex=male)
    sex.a <- det[which(is.finite( det$abdomen ) & det$sex==male),]
    if ( !is.na(sex.a$trip[1])) sex.a$error <- 'sex.a'

    #Sex.c: Indeterminate sex based on measurements taken (chela values where sex=female
    sex.c <- det[which(is.finite( det$chela ) & det$sex==female),]
    if ( !is.na(sex.c$trip[1])) sex.c$error <- 'sex.c'

    #Mat.errors: Unknown Maturity
    mat.errors <- det[which(det$mat ==2 & (is.finite(det$chela) | is.finite(det$abdomen))),]
    if ( !is.na(mat.errors$trip[1])) mat.errors$error <- 'mat.errors'


    names.errors <- list(mat.errors, sex.errors, cw.errors, chela.errors, abdomen.errors, mass.errors, sex.a, sex.c)

    errors = NULL
    for (e in names.errors) {
        if (nrow(e) > 0) errors <- rbind(errors, e[,names(errors)])
    }

    ii = grep(yr.errors, errors$trip)  # added error check  as it causes a break
    if (length(ii) > 0) {
        errors.yearly <- errors[ii,]
        ## JC Jan 2022.. not sure what is going on here
        errors <<- errors  #??
        message("check dataframe 'errors' for the errors")
        if ( !is.na(errors.yearly$trip[1]))  {
            print(errors.yearly)
            write.csv(errors.yearly, file=outfile.errors)
            print("Current Year Morphology Errors saved to file")
            print(outfile.errors)
        }
    }

    write.csv(errors, file=outfile.errors2)
    
    message("Morphology Errors (if any) saved to file: \n")

    message(outfile.errors2)
    
    message("\n\nThe errors found include: \n\n")
    
    cat(errors)

    cat("----\n")
    cat("ERROR CODES are as follows: \n
        Mat.errors: Unknown Maturity \n
        Sex.errors: Unknown Sex \n
        Cw.errors: Carapace Width below 5 or greater than 185 \n
        Chela.errors: Chela less than 1 or greater than 50 \n
        Abdomen.errors:Abdomen less than 1 and greater than 66 \n
        Mass.errors: Mass less than 1 or greater than 1500 \n
        Sex.a: Indeterminate sex based on measurements taken (abdomen values where sex=male) \n
        Sex.c: Indeterminate sex based on measurements taken (chela values where sex=female \n"
    ) 

    return()
}
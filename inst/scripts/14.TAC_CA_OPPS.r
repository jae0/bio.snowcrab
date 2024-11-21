#######################################################################################
# The first statement here needs to be run in order to update the TAC tables
#
# The following statements are run as a double check to the CA calculations 
# excel sheet. This should be done as errors in the sheet are often found. Make 
# sure these results match the results of the excel sheet.
#######################################################################################


## Must obtain TAC for each area. 
## https://inter-j02.dfo-mpo.gc.ca/mqr/quotareports/snowcrab?rptyear=2024&rptnote=false&lang=en
## Replace the year section of the above address with the desired year

## SUBTRACT 250 mt from the CFA 24 TAC as this is taken off the top for Millbrook.  
CA.insertTAC(cfa23.tac=3972.96, cfa24.tac=3122.12, Millbrook.tac = 250, nens.tac= 937.635, xxxx.tac= 20, year=2024)



##  Must update budget annually and add the correct values tot the database
CA.insertVariables(price_per_mt = 2970, budget = 940447, year = 2023)

## Obtain the percentages by individual from current Fisheries Management contact
CA.insertCDDcsv() 
## Once this is completed it may list individuals that are not listed in the individuals table
## This happens when people buy into the fishery or individuals/companies chage names
## New individuals will need to be added manually to the individuals table after figureing out 
## their representatives

## These functions due a final check to ensure the tables are adding up and write results to file
CA.writeIndividualContributions(2023)
CA.writePartnerContributions(2023)
CA.writeAreaContributions(2023)


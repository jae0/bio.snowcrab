
management_areal_units = function( mau="region" ) {

    # define the main managemnent areas (mau) in terms of the internal names   
    # mau is also a variable name where applicable (mostly logbook data)
    
    management_units_all = list(
        region = list(
            cfanorth = "N-ENS",
            cfasouth = "S-ENS",
            cfa4x = "4X"
        ),
        subarea = list(
            cfanorth = "CFA 20-22",
            cfa23 = "CFA 23",
            cfa24 = "CFA 24",
            cfa4x = "CFA 4X"
        )
    )

    management_units_selected = management_units_all[[mau]]  
    n = length(management_units_selected)

    maus = list(
        management_units = management_units_all,
        list = management_units_selected,
        internal = names(management_units_selected),  
        n = n,
        labels = unname( unlist( management_units_selected ) ),
        color_map = c("#E69F00", "#56B4E9",  "#CC79A7" , "#D55E00", "#F0E442", "purple", "gray" )[1:n],
        shapes = c(15, 17, 19, 21, 23, 25, 27)[1:n]

    )

    return(maus)
}

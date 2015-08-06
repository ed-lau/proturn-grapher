calcLabelingSite <- function(string){
        
        N <- 0.0
        carb <- 0
        nitr <- 0
        oxyg <- 0 
        sulf <- 0
        
        # Deal with modifications first
        if(grepl("\\(15\\..+?\\)",string) == T)  {N = N + 0.00; carb = carb + 0 ; nitr = nitr + 0; oxyg = oxyg + 1; sulf = sulf + 0}
        if(grepl("\\(42\\..+?\\)",string) == T)  {N = N + 0.00; carb = carb + 2 ; nitr = nitr + 0; oxyg = oxyg + 1; sulf = sulf + 0}
        if(grepl("\\(0\\.984\\)",string) == T)  {N = N + 0.00; carb = carb + 0 ; nitr = nitr - 1; oxyg = oxyg + 1; sulf = sulf + 0}
        if(grepl("\\(79\\..+?\\)",string) == T)  {N = N + 0.00; carb = carb + 0 ; nitr = nitr + 0; oxyg = oxyg + 4; sulf = sulf + 0}
        if(grepl("\\(114\\..+?\\)",string) == T) {N = N + 4.12; carb = carb + 4 ; nitr = nitr + 3; oxyg = oxyg + 1; sulf = sulf + 0}
        
        string <- gsub( " *\\(.*?\\) *", "", string)
        
        # Remove one oxygen per peptide bond
        oxyg <- oxyg - nchar(string) +1
        
        # Each character - add labeling sites from Commerford et al., add atoms.
        for (i in 1:nchar(string)){
                char <- substring(string,i,i)
                if(char == "A"){N = N + 4.00; carb = carb + 3 ; nitr = nitr + 1; oxyg = oxyg + 2; sulf = sulf + 0}
                if(char == "C"){N = N + 1.62; carb = carb + 3 ; nitr = nitr + 1; oxyg = oxyg + 2; sulf = sulf + 1}
                if(char == "D"){N = N + 1.89; carb = carb + 4 ; nitr = nitr + 1; oxyg = oxyg + 4; sulf = sulf + 0}
                if(char == "E"){N = N + 3.95; carb = carb + 5 ; nitr = nitr + 1; oxyg = oxyg + 4; sulf = sulf + 0}
                if(char == "F"){N = N + 0.32; carb = carb + 9 ; nitr = nitr + 1; oxyg = oxyg + 2; sulf = sulf + 0}
                if(char == "G"){N = N + 2.06; carb = carb + 2 ; nitr = nitr + 1; oxyg = oxyg + 2; sulf = sulf + 0}
                if(char == "H"){N = N + 2.88; carb = carb + 6 ; nitr = nitr + 3; oxyg = oxyg + 2; sulf = sulf + 0}
                if(char == "I"){N = N + 1.00; carb = carb + 6 ; nitr = nitr + 1; oxyg = oxyg + 2; sulf = sulf + 0}
                if(char == "K"){N = N + 0.54; carb = carb + 6 ; nitr = nitr + 2; oxyg = oxyg + 2; sulf = sulf + 0}
                if(char == "L"){N = N + 0.69; carb = carb + 6 ; nitr = nitr + 1; oxyg = oxyg + 2; sulf = sulf + 0}
                if(char == "M"){N = N + 1.12; carb = carb + 5 ; nitr = nitr + 1; oxyg = oxyg + 2; sulf = sulf + 1}
                if(char == "N"){N = N + 1.89; carb = carb + 4 ; nitr = nitr + 2; oxyg = oxyg + 3; sulf = sulf + 0}
                if(char == "P"){N = N + 2.59; carb = carb + 5 ; nitr = nitr + 1; oxyg = oxyg + 2; sulf = sulf + 0}
                if(char == "Q"){N = N + 3.95; carb = carb + 5 ; nitr = nitr + 2; oxyg = oxyg + 3; sulf = sulf + 0}
                if(char == "R"){N = N + 3.34; carb = carb + 6 ; nitr = nitr + 4; oxyg = oxyg + 2; sulf = sulf + 0}
                if(char == "S"){N = N + 2.61; carb = carb + 3 ; nitr = nitr + 1; oxyg = oxyg + 3; sulf = sulf + 0}
                if(char == "T"){N = N + 0.20; carb = carb + 4 ; nitr = nitr + 1; oxyg = oxyg + 3; sulf = sulf + 0}
                if(char == "V"){N = N + 0.56; carb = carb + 5 ; nitr = nitr + 1; oxyg = oxyg + 2; sulf = sulf + 0}
                if(char == "W"){N = N + 0.08; carb = carb + 11; nitr = nitr + 2; oxyg = oxyg + 2; sulf = sulf + 0}
                if(char == "Y"){N = N + 0.42; carb = carb + 9 ; nitr = nitr + 1; oxyg = oxyg + 3; sulf = sulf + 0}
        }
        
        N <- round(N)
        return(c(N,carb,nitr,oxyg,sulf))
}


# Natural abundance of isotopes from NIST and http://www.iupac.org/publications/pac/83/2/0397/
calc_thr_a <- function(string){
        atoms <- calcLabelingSite(string)
        c <- atoms[2]
        n <- atoms[3]
        o <- atoms[4]
        s <- atoms[5]
        a <- 0.9893^c * 0.99636^n * 0.99757^n * 0.9499^s
        
        return(a)
}



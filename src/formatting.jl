function export_csv(df; filename="symCF_output"::String)
    df_clean=cleanLabels(deepcopy(df))
    CSV.write("$filename.csv",df_clean,header=false)
end

"""
        function export_symbolic_format(net, df;
            inheritancecorrelation=0.0::Float64,
            filename="symCF_output"::String,
            csv=true::Bool,
            macaulay2=false::Bool,
            matlab=false::Bool,
            multigraded=false::Bool,
            singular=false::Bool
            )

# Description
- Print symbolic CFs to csv file. (default)
- Format symbolic CF formulas for use in macaulay2, Matlab, and Singular.
"""
function export_symbolic_format(net, df;
    inheritancecorrelation=0.0::Float64,
    filename="symCF_output"::String,
    csv=true::Bool,
    macaulay2=false::Bool,
    matlab=false::Bool,
    multigraded=false::Bool,
    singular=false::Bool
    )
    
    if(csv) export_csv(df; filename=filename) end
    
    #macaulay2 output
    numCFs=size(df)[1]
    dataframe=deepcopy(df)
    
    params=gettingSymbolicInput(net, dataframe, inheritancecorrelation) 
    
    function rationalize(x;sigdigits=16)
        rational=string(Int(round(x*10^(sigdigits-1),digits=0))//10^(sigdigits-1))
        rational=replace(rational, r"//" => "/")
        return "($rational)"
    end
    
    if(macaulay2)
        filename1=filename*"_macaulay2.m2"
        open(filename1, "w") do file
            str="R = QQ["
            for par in params str=str*par*"," end
            str=str*"C_1..C_$numCFs]\n"
            str=str*"I = ideal(\n"
            i=1
            while i<numCFs
                str1=(dataframe[i,2])
                str1=replace(str1,string(inheritancecorrelation)=>rationalize(inheritancecorrelation))
                str=str*"$str1-C_$i,\n"
                i+=1
            end
            str1=(dataframe[numCFs,2])
            str1=replace(str1,string(inheritancecorrelation)=>rationalize(inheritancecorrelation))
            str=str*"$(str1)-C_$numCFs);\n"
            str=str*"G = eliminate (I, {"
            numparams=length(params)   
            for par in 1:(length(params)-1) str=str*params[par]*"," end
            str=str*"$(params[numparams])})\n"
            str=str*"S = QQ[C_1..C_$numCFs]\nJ = sub(G,S)\ndim J"
            write(file, str)
        end
    end
    
    #matlab output
    if(matlab)
        filename1=filename*"_matlab"
        #str*="Write Matlab file: "
        #str*=(matlab ? "on\n" : "off\n")
        open("$filename1.m", "w") do file 
            str="clc\n"
            str=str*"% Declare variables\n"
            str=str*"syms "
            for par in params str=str*par*" " end
            for i in 1:numCFs str=str*"C$i " end
            str=str*"\n\n%%% Got \"Unrecognized function or variable 'syms'\" error msg? Uncomment the next line to fix it\n"
            str=str*"% restoredefaultpath"
            str=str*"\n\n% matrix of generating polynomials\n"
            str=str*"F=["
            i=1
            while i<numCFs
                str1=(dataframe[i,2])
                str1=replace(str1,string(inheritancecorrelation)=>rationalize(inheritancecorrelation))                
                str=str*"$(str1)-C$i,\n"
                i+=1
            end
            str1=(dataframe[numCFs,2])
            str1=replace(str1,string(inheritancecorrelation)=>rationalize(inheritancecorrelation))            
            str=str*"$(str1)-C$numCFs];\n"
#            str=str*"\n% matrix of generating polynomials\n"
            str=str*"\n% Array of all variables\n"
            str=str*"V=["
            for par in params str=str*par*" " end
            for i in 1:numCFs-1 str=str*"C$i " end
            str=str*"C$numCFs]\n"
            str=str*"\n% Compute dimension\nCoalDim(F,V)"
            str=str*"\n
%% The function that numerically computes the dimension of the variety defined by our CF parametriztion
function CoalDim(c, d)
    n_equations = length(c); % number of equations (or number of C's in the ideal)
    n_parameters = length(d) - n_equations;% number of parameters in X's and R's
    a = n_parameters;
    b = n_equations;
    JF = jacobian(c, d);
    B = [zeros(b,a) -eye(b)];
    rk = rank([JF; B]); % nxn that varnish (Check)
    sz = size([JF; B]); % The size function returns the dimensions of an array (we take the second entry)
    \"The dimension of the variety is \" + string(a - (sz(2) - rk))
end"
            write(file, str)
        end
    end
    
    if(multigraded)
        #filename*="_macaulay"
        open("$(filename)_im.m2", "w") do file
            str="needsPackage \"MultigradedImplicitization\"\n"
            str*="R = QQ["
            for par in params str=str*par*"," end
            str*="T]\n"
            str*="S = QQ["
            str=str*"C_0..C_$numCFs]\n"
            str=str*"im = {\n"
            i=1
            while i<numCFs
                str1=(dataframe[i,2])
                str1=replace(str1,string(inheritancecorrelation)=>rationalize(inheritancecorrelation))                 
                str=str*"$(str1),\n"
                i+=1
            end
            str1=(dataframe[numCFs,2])
            str1=replace(str1,string(inheritancecorrelation)=>rationalize(inheritancecorrelation))            
            str=str*"$(str1)};\n"
            
            str=str*"im = {T} | apply(im, f -> f * T);\n"
            str=str*"phi = map(R, S, im)\n"
            str=str*"d=1\n"
            str=str*"I = time componentsOfKernel(d, phi)\n"
            str=str*"L = flatten values I"
            
            write(file, str)
        end
    end
    
    if(singular)
        open("$filename.sing", "w") do file
            str="ring R = 0, ("
            for par in params str=str*par*"," end
            for i in 1:numCFs str=str*"C$i," end
            str=chop(str)
            str=str*"), dp;\n"
            str=str*"ideal I = \n"
            i=1
            while i<numCFs
                str1=(dataframe[i,2])
                str1=replace(str1,string(inheritancecorrelation)=>rationalize(inheritancecorrelation))               
                str=str*"$(str1)-C$i,\n"
                i+=1
            end
            str1=(dataframe[numCFs,2])
            str1=replace(str1,string(inheritancecorrelation)=>rationalize(inheritancecorrelation))  
            str=str*"$(str1)-C$numCFs;\n"
            
            str=str*"ideal G = eliminate (I, "
            for par in params str=str*par*"*" end
            str=chop(str)
            str=str*");\n"
            
            str=str*"ring S = 0, ("
            for i in 1:numCFs str=str*"C$i," end
            str=chop(str)
            str=str*"), dp;\n"
            str=str*"ideal J = imap(R, G);\n"
            str=str*"int dimJ = dim(J);\n"
            str=str*"dimJ;"
            write(file, str)
        end
    end
    
    return nothing
end
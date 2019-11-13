import os
import Bio.PDB
import csv
import pandas as pd
import numpy as np
pd.set_option('display.max_rows', 500)

def itertest(input_df):
    for index,row in input_df.iterrows():
       print(index,row['phi'],row['psi'])

def phitranslate(arr):
    translatearr = []
    index = 0
    for item in arr:
        if (item < -30.0 and item > -150.0): 
            translatearr.append(-1) 
            index+=1
        elif (item > 30.0 and item < 150.0):    
            translatearr.append(1)
            index+=1
        elif index == 0: #what if the -30 to 30 is the first value? - just make it -1? or make it the same as the next value to kill alternation
            if (arr[1] < -30.0 and arr[1] > -150.0):
                translatearr.append(-1)
                index+=1
            elif(arr[1] > 30.0 and arr[1] < 150.0):
                translatearr.append(1)
                index+=1
        else: #what to do if between -30 and 30? - should probably make those the same as previous value to kill alternation
            translatearr.append(translatearr[-1])
    return translatearr

def phitranslate2(arr):
    translatearr = []
    index = 0
    for item in arr:
        if (item < -30.0 and item > -150.0): 
            translatearr.append(-1) 
            index+=1
        elif (item > 30.0 and item < 150.0):    
            translatearr.append(1)
            index+=1
        elif index == 0: #what if the -30 to 30 is the first value? - just make it -1? or make it the same as the next value to kill alternation
            if (arr[1] < -30.0 and arr[1] > -150.0):
                translatearr.append(-1)
                index+=1
            elif(arr[1] > 30.0 and arr[1] < 150.0):
                translatearr.append(1)
                index+=1
            else:
                translatearr.append(0)
        else: #what to do if between -30 and 30? - should probably make those the same as previous value to kill alternation
            #translatearr.append(translatearr[-1])
            translatearr.append(0)
            index+=1
    return translatearr



def phipsitranslate2(input_df):
    translatearr = []
    for index,row in input_df.iterrows():
        try:
            if ((row['phi'] < -30.0 and row['phi'] > -150.0) and (row['psi'] < 45.0 and row['psi'] > -90.0)): #right
                translatearr.append(-1)             
            elif ((row['phi'] < 150.0 and row['phi'] > 30.0) and (row['psi'] < 90.0 and row['psi'] > -45.0)): #left
                translatearr.append(1)
            elif (index == 0): #what if non-alpha is the first value? or make it the same as the next value to kill alternation
                if (input_df.iloc[1]['phi'] < -30.0 and input_df.iloc[1]['phi'] > -150.0) and (input_df.iloc[1]['psi'] < 45.0 and input_df.iloc[1]['phi'] > -90.0): #right
                    translatearr.append(-1)
              
                elif(input_df.iloc[1]['phi'] < 150.0 and input_df.iloc[1]['phi'] > 30.0) and (input_df.iloc[1]['psi'] < 90.0 and input_df.iloc[1]['psi'] > -45.0): #left
                    translatearr.append(1)
                else:
                    translatearr.append(0)
              
            else: #what if non-alpha? - should probably make those the same as previous value to kill alternation
                #translatearr.append(translatearr[-1]) #causing errors, need to make non-alpha values 0 and modify zigzags
                translatearr.append(0)
        except IndexError:
            continue
    return translatearr 

def zigzags(input):
    #input = iter(phipsitranslate(input))
    input = iter(phitranslate(input['phi']))
    #input = iter(np.sign(df['phi']))
    stack = None
    index = 0
    try:
        stack = [[index,next(input)]]
        index+=1
        while True:
            if len(stack) < 2:
                stack.append([index,next(input)])
                index+=1
            else:
                stack = stack[-2:]
            a, b = stack
            if(a[1] == b[1]):
                #yield (a,)
                stack = [b]
                continue
            zig = a[1] > b[1]
            while True:
                prev = stack[-1]
                this = [index,next(input)]
                index+=1
                if prev[1] == this[1] or zig == (prev[1] > this[1]):
                    break
                stack.append(this)
                zig = not zig
            yield tuple(stack)
            stack.append(this)
    except StopIteration:
        pass
    if stack:
        yield tuple(stack)


def zigzags2(input):
    #input = iter(phipsitranslate2(input))
    input = iter(phitranslate2(input['phi']))
    ##input = iter(np.sign(df['phi']))
    stack = None
    index = 0
    try:
        stack = [[index,next(input)]]
        index+=1
        while True:
            if len(stack) < 2:
                stack.append([index,next(input)])
                index+=1
            else:
                stack = stack[-2:]
            a, b = stack
            if((a[1] == b[1]) or (a[1] ==0) or (b[1]==0)):
                yield (a,)
                stack = [b]
                continue
            zig = a[1] > b[1]
            while True:
                prev = stack[-1]
                this = [index,next(input)]
                index+=1
                if prev[1] == this[1] or zig == (prev[1] > this[1]) or  (prev[1] ==0) or (this[1]==0) :
                    break
                stack.append(this)
                zig = not zig
            yield tuple(stack)
            stack.append(this)
    except StopIteration:
        pass
    if stack:
        yield tuple(stack)


def zigzag(df_original,input,df_results):
    bigzigs = list(input)
    for lis in bigzigs:
        if len(lis) > 3:
            startres = lis[0][0]
            #print(df.loc[startres,'name_chain_res'])
            #print(len(lis))
            #print(lis)
            longlist = [startres+1,df_original.loc[startres]['name_chain_res'],len(lis)]
            df_results.loc[len(df_results)] = longlist
    return df_results


def runzig():
    df_results = pd.DataFrame(columns = ['line_number','name_chain_res','length'])
    os.chdir("C:/Users/Shaheer Rizwan/Documents/ramachandran/phipsi_xray_chains")
    count = 0
    for filename in os.listdir("C:/Users/Shaheer Rizwan/Documents/ramachandran/phipsi_xray_chains"):
        with open(filename,"r") as tsvfile:
            df = pd.read_csv(filename, sep='\t',usecols=[0,1,2], names=['name_chain_res','phi','psi'])
            tsvfile.close()
            print("%s - %s out of 14,998 files" % (filename,count))
            count+=1
            #print(phipsitranslate(df))
            #print(df)
            bigzig = zigzags2(df)
            df_results = zigzag(df,bigzig,df_results)
        
    os.chdir("C:/Users/Shaheer Rizwan/Documents/ramachandran/")
    df_results.to_csv("test_results_phipsi_xray_chains_rev4.csv")
    print("DONE")

#df_results = pd.DataFrame(columns = ['line_number','name_chain_res','length'])
#os.chdir("C:/Users/Shaheer Rizwan/Documents/ramachandran/phipsi_xray_chains")
#df = pd.read_csv('pdb1b35_chainB_biopython.tsv', sep='\t',usecols=[0,1,2], names=['name_chain_res','phi','psi'])
#df['alpha'] = phitranslate2(df['phi'])
#print(df)
#bigzig = zigzags2(df)
#print(zigzag(df,bigzig,df_results))

runzig()
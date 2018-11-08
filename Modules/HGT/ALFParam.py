'''Generate ALFSIM param files to simulate evolution of HGT genes'''


def get_options():
    import argparse

    description = 'Generate ALFSIM param files to simulate evolution of HGT genes'
    parser = argparse.ArgumentParser(description=description)

    parser.add_argument('--ALFParams',
                        help='the ALFParameter file from which the sub/indel/rateHetr parameters to be imported')

    parser.add_argument('--HGTRate',
                        default='10',
                        help="Frequency of Horizontal Gene Transfer events [Default: 10]")
    parser.add_argument('--Output',
                        help="Output file including paramters to be used by ALFSIM for simulation along with destination of ALF simulation data")


    return (parser.parse_args())


options = get_options()


#configfile: "/home/masih/Baterial_simulator/Pipeline/ConfigFile.yaml"
import os



def ALFParam (ALFParams,HGTRate,Output):
    ParamDict={}
    with open (ALFParams,'r') as file:
        line=file.readline().strip()
        while line:
            if ':=' in line:
                ParamDict[line.split(':=')[0].strip()]=line.split(':=')[1]
            line=file.readline().strip()
    print(ParamDict)            
    for HGTs in range(int(HGTRate)):
        txt=open(os.path.join(Output,str(HGTs)+'.drw'),'w')
        txt.write('SetRand(2345):\n')
        txt.write('%s := %s\n' % ('realorganism',"'"+os.path.join(Output,str(HGTs)+'.db')+"' ;"))
        txt.write('%s := %s\n' % ('treeType',ParamDict['treeType']))
        txt.write('%s := %s\n' % ('treeFile',ParamDict['treeFile']))
        txt.write('%s := %s\n' % ('unitIsPam',ParamDict['unitIsPam']))
        txt.write('%s := %s\n' % ('substModels',ParamDict['substModels']))
        txt.write('%s := %s\n' % ('IndelModels',ParamDict['IndelModels']))
        txt.write('%s := %s\n' % ('rateVarModels',ParamDict['rateVarModels']))
        txt.write('%s := %s\n' % ('mname',ParamDict['mname']))
        txt.write('%s := %s\n' % ('wdir',"'"+os.path.join(Output,str(HGTs)+'/')+"' ;"))
        txt.write("simOutput := {\n'Newick',\n'Fasta',\nNULL}:")
        txt.close()

ALFParam(options.ALFParams,options.HGTRate,options.Output)
        
        
        
    
    
        
        
from dflow import (
    Executor,
    ShellOPTemplate,
    InputArtifact,
    InputParameter,
    OutputArtifact,
    OutputParameter,
    Workflow,
    Step,
    Steps,
    argo_range,
    argo_sequence,
    upload_artifact,
    download_artifact
)
from dflow.python import (
    PythonOPTemplate,
    OP,
    OPIO,
    OPIOSign,
    Artifact,
    Slices,
    upload_packages
)
from numpy import array, double, float64
if "__file__" in locals():
    upload_packages.append(__file__)
from pathlib import Path
import time
from typing import List
import urllib3
urllib3.disable_warnings(urllib3.exceptions.InsecureRequestWarning)
import ast
import numpy as np



class ReadCommonDict(OP):
    def __init__(self):
        pass

    @classmethod
    def get_input_sign(cls):
        return OPIOSign({
            'commonDict_path' : Artifact(Path),
            'commonDict_py' : Artifact(Path)
        })

    @classmethod
    def get_output_sign(cls):
        return OPIOSign({
            'cbDict' : dict(work_dir=str,output_fln=str,chemMech=str,
            transModel=str,scaled_PV=bool,f_min=float,f_max=float,
            nchemfile=int,nVarCant=int,CASENAME=str,p=float,Lx=float,
            fuel_species=str,fuel_C=float,fuel_H=float,T_fuel=float,
            X_fuel=float,T_ox=float,X_ox=str,nscal_BC=int,cUnifPts=int,
            n_points_z=int,n_points_c=int,n_points_h=int,nYis=int,
            spc_names=list,int_pts_z=int,int_pts_c=int,int_pts_gz=int,
            int_pts_gc=int,int_pts_gcz=int,nScalars=int,small=float,
            n_procs=int,solIdx=int,nSpeMech=int,z_space=float64,
            c_space=float64,z=float64,c=float64,gz=float64,gc=float64,
            gcz=float64),
        })

    @OP.exec_sign_check
    def execute(
            self,
            op_in : OPIO,
    ) -> OPIO:
        from os import path,mkdir
        pyPath=path.dirname(op_in["commonDict_py"])
        import sys
        sys.path.append(pyPath)
        import ast
        import read_commonDict
        import numpy as np

        cbDict=read_commonDict.read_commonDict(op_in["commonDict_path"])

        #-----done read cbDict--------#


        return OPIO({
            "cbDict": cbDict
        })



class CanteraSim(OP):
    def __init__(self):
        pass

    @classmethod
    def get_input_sign(cls):
        return OPIOSign({
            'cbDict' : dict(work_dir=str,output_fln=str,chemMech=str,
            transModel=str,scaled_PV=bool,f_min=float,f_max=float,
            nchemfile=int,nVarCant=int,CASENAME=str,p=float,Lx=float,
            fuel_species=str,fuel_C=float,fuel_H=float,T_fuel=float,
            X_fuel=float,T_ox=float,X_ox=str,nscal_BC=int,cUnifPts=int,
            n_points_z=int,n_points_c=int,n_points_h=int,nYis=int,
            spc_names=list,int_pts_z=int,int_pts_c=int,int_pts_gz=int,
            int_pts_gc=int,int_pts_gcz=int,nScalars=int,small=float,
            n_procs=int,solIdx=int,nSpeMech=int,z_space=float64,
            c_space=float64,z=float64,c=float64,gz=float64,gc=float64,
            gcz=float64),
            'canteraSim_py' : Artifact(Path),
            # "phi_i" : int
        })

    @classmethod
    def get_output_sign(cls):
        return OPIOSign({
            'canteraSimPath' : Artifact(str)
        })

    @OP.exec_sign_check
    def execute(
            self,
            op_in : OPIO,
    ) -> OPIO:
        from os import path
        pyPath=path.dirname(op_in["canteraSim_py"])
        import sys
        sys.path.append(pyPath)
        import multiProcessing_canteraSim

        return OPIO({
            "canteraSimPath": multiProcessing_canteraSim.canteraSim(op_in["cbDict"])
        })


class InterpToMeshgrid(OP):
    def __init__(self):
        pass

    @classmethod
    def get_input_sign(cls):
        return OPIOSign({
            'cbDict' : dict(work_dir=str,output_fln=str,chemMech=str,
            transModel=str,scaled_PV=bool,f_min=float,f_max=float,
            nchemfile=int,nVarCant=int,CASENAME=str,p=float,Lx=float,
            fuel_species=str,fuel_C=float,fuel_H=float,T_fuel=float,
            X_fuel=float,T_ox=float,X_ox=str,nscal_BC=int,cUnifPts=int,
            n_points_z=int,n_points_c=int,n_points_h=int,nYis=int,
            spc_names=list,int_pts_z=int,int_pts_c=int,int_pts_gz=int,
            int_pts_gc=int,int_pts_gcz=int,nScalars=int,small=float,
            n_procs=int,solIdx=int,nSpeMech=int,z_space=float64,
            c_space=float64,z=float64,c=float64,gz=float64,gc=float64,
            gcz=float64),
            'interp_py' : Artifact(Path),
            'canteraSimPath': Artifact(str),
        })

    @classmethod
    def get_output_sign(cls):
        return OPIOSign({
            'interpPath' : Artifact(str)
        })

    @OP.exec_sign_check
    def execute(
            self,
            op_in : OPIO,
    ) -> OPIO:
        from os import path
        pyPath=path.dirname(op_in["interp_py"])
        import sys
        sys.path.append(pyPath)
        import interpToMeshgrid


        return OPIO({
            "interpPath": interpToMeshgrid.interpLamFlame(op_in["cbDict"],
                    op_in["canteraSimPath"])
        })


class PDFSim(OP):
    def __init__(self):
        pass

    @classmethod
    def get_input_sign(cls):
        return OPIOSign({
            'cbDict' : dict(work_dir=str,output_fln=str,chemMech=str,
            transModel=str,scaled_PV=bool,f_min=float,f_max=float,
            nchemfile=int,nVarCant=int,CASENAME=str,p=float,Lx=float,
            fuel_species=str,fuel_C=float,fuel_H=float,T_fuel=float,
            X_fuel=float,T_ox=float,X_ox=str,nscal_BC=int,cUnifPts=int,
            n_points_z=int,n_points_c=int,n_points_h=int,nYis=int,
            spc_names=list,int_pts_z=int,int_pts_c=int,int_pts_gz=int,
            int_pts_gc=int,int_pts_gcz=int,nScalars=int,small=float,
            n_procs=int,solIdx=int,nSpeMech=int,z_space=float64,
            c_space=float64,z=float64,c=float64,gz=float64,gc=float64,
            gcz=float64),
            'pdfSim_py' : Artifact(Path),
            'interpPath': Artifact(str),
        })

    @classmethod
    def get_output_sign(cls):
        return OPIOSign({
            'pdfPath' : Artifact(str)
        })

    @OP.exec_sign_check
    def execute(
            self,
            op_in : OPIO,
    ) -> OPIO:
        from os import path
        pyPath=path.dirname(op_in["pdfSim_py"])
        import sys
        sys.path.append(pyPath)
        import PDF_Sim_multiProcessing


        return OPIO({
            "pdfPath": PDF_Sim_multiProcessing.PDF_Intergrate(op_in["cbDict"],
                    op_in["interpPath"])
        })


class Assemble(OP):
    def __init__(self):
        pass

    @classmethod
    def get_input_sign(cls):
        return OPIOSign({
            'cbDict' : dict(work_dir=str,output_fln=str,chemMech=str,
            transModel=str,scaled_PV=bool,f_min=float,f_max=float,
            nchemfile=int,nVarCant=int,CASENAME=str,p=float,Lx=float,
            fuel_species=str,fuel_C=float,fuel_H=float,T_fuel=float,
            X_fuel=float,T_ox=float,X_ox=str,nscal_BC=int,cUnifPts=int,
            n_points_z=int,n_points_c=int,n_points_h=int,nYis=int,
            spc_names=list,int_pts_z=int,int_pts_c=int,int_pts_gz=int,
            int_pts_gc=int,int_pts_gcz=int,nScalars=int,small=float,
            n_procs=int,solIdx=int,nSpeMech=int,z_space=float64,
            c_space=float64,z=float64,c=float64,gz=float64,gc=float64,
            gcz=float64),
            'assemble_py' : Artifact(Path),
            'pdfPath': Artifact(str),
        })

    @classmethod
    def get_output_sign(cls):
        return OPIOSign({
            'flareTabPath' : Artifact(str)
        })

    @OP.exec_sign_check
    def execute(
            self,
            op_in : OPIO,
    ) -> OPIO:
        from os import path
        pyPath=path.dirname(op_in["assemble_py"])
        import sys
        sys.path.append(pyPath)
        import assemble


        return OPIO({
            "flareTabPath": assemble.assemble(op_in["cbDict"],
                    op_in["pdfPath"])
        })

#serial
class importCase(OP):
    def __init__(self):
        pass

    @classmethod
    def get_input_sign(cls):
        return OPIOSign({
            'case_path' : Artifact(str),
            'flareTab_path' : Artifact(str),
            'Ubulk' : float,
            'caseIndex' : int,
            'n_proc_mpirun' : int,
        })

    @classmethod
    def get_output_sign(cls):
        return OPIOSign({
            'zip_path' : Artifact(str),
        })

    @OP.exec_sign_check
    def execute(
            self,
            op_in : OPIO,
    ) -> OPIO:
        import zipfile
        import os,shutil
        import fileinput
        import sys

        # modify case input
        casePath = op_in["case_path"]
        replaced_file = casePath + "/0/U"

        def replacement(file, previousw, nextw):
            for line in fileinput.input(file, inplace=1):
                line = line.replace(previousw, nextw)
                sys.stdout.write(line)
        
        var_old = "49.6"
        var_new = str(op_in["Ubulk"])
        replacement(replaced_file,var_old,var_new)

        # modify decomposePar mesh block Number
        replaced_file = casePath + "/system/decomposeParDict"
        var_old = "numberOfSubdomains 4;"
        var_new = "numberOfSubdomains " + str(op_in["n_proc_mpirun"]) + ";"
        replacement(replaced_file,var_old,var_new)

        shutil.move((op_in["flareTab_path"]),casePath)

        zip_path="/tmp/case_"+str(op_in["caseIndex"])+".zip"

        with zipfile.ZipFile(zip_path,'w') as target:
            for i in os.walk(casePath):
                for n in i[2]:
                    target.write(''.join((i[0],"/",n)))
        
        print(os.listdir("/tmp"))

        return OPIO({
            'zip_path': zip_path
        })

#parallel
class runCase(OP):
    def __init__(self):
        pass

    @classmethod
    def get_input_sign(cls):
        return OPIOSign({
            'case_path' : Artifact(str),
            'flareTab_path' : Artifact(str),
            'Ubulk' : float,
            'n_proc_mpirun' : int,
            'caseIndex' :str,
        })

    @classmethod
    def get_output_sign(cls):
        return OPIOSign({
            'zip_path' : Artifact(str),
        })

    @OP.exec_sign_check
    def execute(
            self,
            op_in : OPIO,
    ) -> OPIO:
        import zipfile
        import os,shutil
        import fileinput
        import sys
        import subprocess

        # modify case input
        casePath = op_in["case_path"]
        replaced_file = casePath + "/0/U"

        def replacement(file, previousw, nextw):
            for line in fileinput.input(file, inplace=1):
                line = line.replace(previousw, nextw)
                sys.stdout.write(line)
        
        var_old = "49.6"
        var_new = str(op_in["Ubulk"])
        replacement(replaced_file,var_old,var_new)

        # modify decomposePar mesh block Number
        replaced_file = casePath + "/system/decomposeParDict"
        var_old = "numberOfSubdomains 4;"
        var_new = "numberOfSubdomains " + str(op_in["n_proc_mpirun"]) + ";"
        replacement(replaced_file,var_old,var_new)

        shutil.move((op_in["flareTab_path"]),casePath)
        subcasePath = "./" + str(op_in["caseIndex"] + "/")
        shutil.copytree(casePath,subcasePath)

        subprocess.call(["bash","-c","cd " + str(subcasePath) + \
                "&& source $HOME/OpenFOAM/OpenFOAM-7/etc/bashrc && blockMesh \
                && decomposePar && mpirun -np " + str(op_in["n_proc_mpirun"]) + \
                    " --allow-run-as-root --oversubscribe RANSflareFoam -parallel"])

        zip_path = str(op_in["caseIndex"]) +".zip"

        with zipfile.ZipFile(zip_path,'a') as target:
            for i in os.walk(subcasePath):
                for n in i[2]:
                    target.write(''.join((i[0],"/",n)))
        
        return OPIO({
            'zip_path': zip_path
        })


if __name__ == "__main__":  
    start=time.time()

    #---------------- WorkFlow ----------------#

    #************** Case SetUp **************#
    # Main jet velocity of Sandia Flame C/D/E/F
    # Ubulk = [29.7,49.6,74.4,99.2]
    # caseIndex = ["Sandia_Flame_C","Sandia_Flame_D","Sandia_Flame_E","Sandia_Flame_F"]
    Ubulk = [29.7,49.6]  # to-be-changed parameter list
    caseIndex = ["Sandia_Flame_C","Sandia_Flame_D"]  #case names
    n_proc_mpirun = 7  #processors per case
    artifact_case=upload_artifact("SandiaFlame_RANS/")  #upload case folder
    parallel_cases=True   #if run CFD cases simultaneously
    #************* End Case SetUp ************#

    wf = Workflow(name="flare")

    artifact_commonDict0=upload_artifact("commonDict.txt")
    artifact_commonDict1=upload_artifact("read_commonDict.py")
    load_dict = Step(name="loadCommonDict",template=PythonOPTemplate(ReadCommonDict,
                        image="wangdongchina/pythoncantera:v0.4",command=["ipython"],
                        image_pull_policy="IfNotPresent",OutputParameter={"cbDict"}),
                    artifacts={"commonDict_path": artifact_commonDict0,
                        "commonDict_py":artifact_commonDict1},key="load_dict")
    wf.add(load_dict)


    artifact_canteraSim=upload_artifact("multiProcessing_canteraSim.py")
    cantera_Sim=Step(name='canteraSim',template=PythonOPTemplate(CanteraSim,
                        image="wangdongchina/pythoncantera:v0.4",command=["ipython"],
                        image_pull_policy="IfNotPresent",OutputArtifact="canteraSimPath"),
                    parameters={"cbDict":load_dict.outputs.parameters["cbDict"]},
                    artifacts={"canteraSim_py":artifact_canteraSim},key="cantera_Sim")
    wf.add(cantera_Sim)

    artifact_interpToMesh=upload_artifact("interpToMeshgrid.py")
    interp_Meshgrid=Step("interpToMeshgrid",template=PythonOPTemplate(InterpToMeshgrid,
                            image="wangdongchina/pythoncantera:v0.4",command=["ipython"],
                            image_pull_policy="IfNotPresent",OutputArtifact="interpPath"),
                        parameters={"cbDict":load_dict.outputs.parameters["cbDict"]},
                        artifacts={"interp_py":artifact_interpToMesh,
                            "canteraSimPath":cantera_Sim.outputs.artifacts["canteraSimPath"]},
                        key="interp_Meshgrid")
    wf.add(interp_Meshgrid)

    artifact_pdfSim = upload_artifact("PDF_Sim_multiProcessing.py")
    pdf_Sim = Step(name="pdfSim", template=PythonOPTemplate( PDFSim,
                        image="wangdongchina/pythoncantera:v0.4",command=["ipython"],
                        image_pull_policy="IfNotPresent",OutputArtifact="pdfPath"),
                    parameters={"cbDict":load_dict.outputs.parameters["cbDict"]},
                    artifacts={"pdfSim_py":artifact_pdfSim,
                            "interpPath":interp_Meshgrid.outputs.artifacts["interpPath"]},
                    key="pdf_Sim")
    wf.add(pdf_Sim)

    artifact_assemble = upload_artifact("assemble.py")
    assemble_Flare = Step(name="assembleFlare", template=PythonOPTemplate( Assemble,
                            image="wangdongchina/pythoncantera:v0.4",command=["ipython"],
                            image_pull_policy="IfNotPresent",OutputArtifact="flareTabPath"),
                        parameters={"cbDict":load_dict.outputs.parameters["cbDict"]},
                        artifacts={"assemble_py":artifact_assemble,
                                "pdfPath":pdf_Sim.outputs.artifacts["pdfPath"]},
                        key="assemble_Flare")
    wf.add(assemble_Flare)


    ############### parallel computation for cases ###################
    runcase=PythonOPTemplate(runCase,
            image="wangdongchina/pythoncantera:v0.4",command=["ipython"],
            image_pull_policy="IfNotPresent",OutputArtifact='zip_path',
            slices=Slices("{{item}}",input_parameter=["Ubulk","caseIndex"],
                    output_artifact=["zip_path"]) )
    runCase=Step(name="runcase",template=runcase,when="%s == True" % (parallel_cases),
                artifacts={"case_path": artifact_case,
                    "flareTab_path":assemble_Flare.outputs.artifacts["flareTabPath"]},
                parameters={"Ubulk":Ubulk,"n_proc_mpirun":n_proc_mpirun,"caseIndex":caseIndex},
                with_sequence=argo_sequence(len(Ubulk)),key="runcase")
    wf.add(runCase)

    ############### serial computation for cases ###################
    for caseIndex in range(len(Ubulk)):
        artifact_flareTab=upload_artifact("flare.tbl")
        artifact_case=upload_artifact("SandiaFlame_RANS/")
        importcase=PythonOPTemplate(importCase,
                image="wangdongchina/pythoncantera:v0.4",command=["ipython"],
                image_pull_policy="IfNotPresent",OutputArtifact='zip_path')

        blockmesh=ShellOPTemplate(
            image="wangdongchina/pythoncantera:v0.4",
            command = ["bash"],
            script="ls -R /tmp && cd /tmp/case_{{inputs.parameters.caseIndex}}.tar/tmp && ls -R \
                && unzip case_{{inputs.parameters.caseIndex}}.zip \
                && cd ./tmp/inputs/artifacts/case_path/SandiaFlame_RANS \
                && source $HOME/OpenFOAM/OpenFOAM-7/etc/bashrc \
                && blockMesh && cd .. \
                && tar -cvf blockMesh_{{inputs.parameters.caseIndex}}.tar SandiaFlame_RANS \
                && mv ./blockMesh_{{inputs.parameters.caseIndex}}.tar /tmp"
            )
        blockmesh.inputs.parameters={"caseIndex":InputParameter()}
        blockmesh.inputs.artifacts= \
            {"SandiaFlame_RANS":InputArtifact(path='/tmp/case_'+str(caseIndex)+'.tar')}
        blockmesh.outputs.artifacts={"blockMeshFile": 
            OutputArtifact(path="/tmp/blockMesh_"+str(caseIndex)+".tar")}

        decomposepar=ShellOPTemplate(
            image="wangdongchina/pythoncantera:v0.4",
            command = ["bash"],
            script="cd /tmp && tar -xvf blockMesh_{{inputs.parameters.caseIndex}}.tar && cd SandiaFlame_RANS/ \
                && source $HOME/OpenFOAM/OpenFOAM-7/etc/bashrc \
                && decomposePar && cd .. && tar -cvf decomposePar_{{inputs.parameters.caseIndex}}.tar SandiaFlame_RANS"
            )
        decomposepar.inputs.parameters={"caseIndex":InputParameter()}
        decomposepar.inputs.artifacts={"blockMeshFile":InputArtifact(path="/tmp/blockMesh_"+str(caseIndex)+".tar")}
        decomposepar.outputs.artifacts={"decomposeParFile":OutputArtifact(path="/tmp/decomposePar_"+str(caseIndex)+".tar")}

        mpirun=ShellOPTemplate(
            image="wangdongchina/pythoncantera:v0.4",
            command = ["bash"],requests={"cpu": 7},
            script="cd /tmp && tar -xvf decomposePar_{{inputs.parameters.caseIndex}}.tar && cd SandiaFlame_RANS/ \
                && source $HOME/OpenFOAM/OpenFOAM-7/etc/bashrc \
                && mpirun -np {{inputs.parameters.n_proc_mpirun}} --allow-run-as-root --oversubscribe RANSflareFoam -parallel \
                && cd .. && tar -cvf mpirun_{{inputs.parameters.caseIndex}}.tar SandiaFlame_RANS"
            )
        mpirun.inputs.parameters={"caseIndex":InputParameter(),"n_proc_mpirun":InputParameter()}
        mpirun.inputs.artifacts={"decomposeParFile":InputArtifact(path="/tmp/decomposePar_"+str(caseIndex)+".tar")}
        mpirun.outputs.artifacts={"mpirunFile":OutputArtifact(path="/tmp/mpirun_"+str(caseIndex)+".tar")}
    

        ImportCase=Step(name="loadcase"+str(caseIndex),template=importcase,when="%s != True" % (parallel_cases),
                        artifacts={"case_path": artifact_case,
                            "flareTab_path":artifact_flareTab},
                        parameters={"Ubulk":Ubulk[caseIndex],"caseIndex":caseIndex,"n_proc_mpirun":n_proc_mpirun})
        wf.add(ImportCase)
        BlockMesh=Step(name="blockmesh"+str(caseIndex),template=blockmesh,when="%s != True" % (parallel_cases),
                        parameters={"caseIndex":str(caseIndex)},
                    artifacts={"SandiaFlame_RANS":ImportCase.outputs.artifacts["zip_path"]})
        wf.add(BlockMesh)
        DecomposePar=Step(name="decomposepar"+str(caseIndex),template=decomposepar,parameters={"caseIndex":str(caseIndex)},
            artifacts={"blockMeshFile":BlockMesh.outputs.artifacts["blockMeshFile"]},when="%s != True" % (parallel_cases))
        wf.add(DecomposePar)
        MpiRun=Step(name="mpirun"+str(caseIndex),template=mpirun,parameters={"caseIndex":str(caseIndex),"n_proc_mpirun":n_proc_mpirun},
            artifacts={"decomposeParFile":DecomposePar.outputs.artifacts["decomposeParFile"]},when="%s != True" % (parallel_cases))
        wf.add(MpiRun)
    

    wf.submit()

    while wf.query_status() in ["Pending", "Running"]:
        time.sleep(1)

    assert(wf.query_status() == "Succeeded")

    # # ---------reuse step--------#
    # load_dict1 = wf.query_step(key="load_dict")[0]
    # cantera_Sim1 = wf.query_step(key="cantera_Sim")[0]
    # interp_Meshgrid1 = wf.query_step(key="interp_Meshgrid")[0]
    # pdf_Sim1 = wf.query_step(key="pdf_Sim")[0]
    # assemble_Flare1 = wf.query_step(key="assemble_Flare")[0]
    # runcase = wf.query_step(key="runcase")[0]

    # wf = Workflow("reused-flare")
    # wf.submit(reuse_step=[load_dict1,cantera_Sim1,interp_Meshgrid1])
    # #----------------------------#


    end=time.time()
    print('Running time: %s Seconds'%(end-start))

    if (parallel_cases):
        for caseIndex in range(len(Ubulk)):
            queryStep=wf.query_step(name="runcase")[0]
            assert(queryStep.phase == "Succeeded")
            download_artifact(queryStep.outputs.artifacts["zip_path"])
    else:
        for caseIndex in range(len(Ubulk)):
            queryStep=wf.query_step(name="mpirun"+str(caseIndex))[0]
            assert(queryStep.phase == "Succeeded")
            download_artifact(queryStep.outputs.artifacts["mpirunFile"])





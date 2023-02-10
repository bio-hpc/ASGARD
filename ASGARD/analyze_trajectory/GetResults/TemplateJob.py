import os
class TemplateJob():

    def __init__(self, cfg):
        self.cfg = cfg
        self.cores = "5"
        self.time = "12:00:00"
        self.name_job = "g_"+self.cfg.sufijo

        #self.maxWarning = "5"

        self.lst_header = []
        self.get_headers()

    def get_headers(self):
        cmd = "whereis qsub |grep bin"
        gst = self.cfg.tools.execute.run_with_fail(cmd)
        if gst != "failed":
            self.queue_manager="Sun grid"
            proyecto = "nn2855k"
            self.lst_header.append("#!/bin/sh")
            self.lst_header.append("#PBS -A " + proyecto)
            self.lst_header.append("#PBS -o " + self.cfg.folder + "jobs_out/" )
            self.lst_header.append("#PBS -e " + self.cfg.folder + "jobs_out/")
            self.lst_header.append("#PBS -N " + self.name_job)
            self.lst_header.append("#PBS -l walltime=" + self.time)
            self.lst_header.append("#PBS -l nodes=1:ppn=" + self.cores)

            self.run_job = "qsub "
            self.dependency_job_cmd = "qsub -W depend=afterok"
            self.out_job = "##PBS -o "
            self.err_job = "#PBS -e "

        elif (self.cfg.tools.execute.run_with_fail("whereis sbatch |grep bin") != "faild"):
            self.queue_manager = "Slurm"
            self.lst_header.append("#!/bin/sh")
            self.lst_header.append("#SBATCH --output=" + self.cfg.folder + "jobs_out/")
            self.lst_header.append("#SBATCH --error=" + self.cfg.folder + "jobs_out/")
            self.lst_header.append("#SBATCH -J " + self.name_job)
            self.lst_header.append("#SBATCH --time=" + self.time)
            self.lst_header.append("#SBATCH --cpus=" + self.cores)
            self.run_job = "sbatch "
            self.dependency_job_cmd = "sbatch --depend=afterany"
            self.out_job = "#SBATCH --output="
            self.err_job = "#SBATCH --error="

        else:
            print("ERROR: Queue Manager no detected")
            exit()

    def replace_out_path(self, template_file, file):
        """
            remplaza la salida por una personlaizada para ese script
        """
        suf = os.path.splitext(os.path.basename(template_file))[0];
        for i in self.lst_header:
            if self.out_job in i:
                i += suf + ".out"
            if self.err_job in i:
                i += suf + ".err"
            file.write(i + "\n")
        file.write("\n")
        file.write('export GMX_MAXBACKUP=-1\n\n')

    def execute_job(self, script_template_file, lst_cmd):
        """
            Ejecuta un script en paralelo o sucuencial
        """
        file = open(script_template_file, "w")
        self.replace_out_path(script_template_file, file)

        for i in lst_cmd:
            file.write(i + "\n\n")
        if self.cfg.p_sequential:
            cmd = "sh {}".format(script_template_file)
        else:
            cmd = '{} {}'.format(self.run_job, script_template_file)
        file.close()
        out = self.cfg.tools.execute.run(cmd)
        if not self.cfg.p_sequential:  # job dependencies
            self.cfg.lstJobs.append(out.split(" ")[3])

    def execute_job_with_dependency(self, script_template_file, lst_cmd, dependencias):
        file = open(script_template_file, "w")
        self.replace_out_path(script_template_file, file)
        for i in lst_cmd:
            file.write("\n"+i+"\n")
        if self.cfg.p_sequential:
            cmd = "sh "+script_template_file
        else:
            cmd = self.cfg.dependency_job_cmd +":"+dependencias+" "+ script_template_file
        file.close()
        print(self.cfg.tools.execute.run(cmd))
#
# Copyright (C) 2014 INRA
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#

__author__ = 'Frederic Escudie - Maria Bernard'
__copyright__ = 'Copyright (C) 2015 INRA'
__license__ = 'GNU General Public License'
__version__ = '0.2.0'

import os
import sys
import time
import subprocess
from subprocess import Popen, PIPE
import gzip

##################################################################################################################################################
#
# FUNCTIONS
#
##################################################################################################################################################


def is_gzip(file):
    """
    @return: [bool] True if the file is gziped.
    @param file : [str] Path to processed file.
    """
    is_gzip = None
    FH_input = gzip.open( file )
    try:
        FH_input.readline()
        is_gzip = True
    except:
        is_gzip = False
    finally:
        FH_input.close()
    return is_gzip


def which(exec_name):
    """
    @summary: Returns the software absolute path.
    @param exec_name: [str] The software file (example : blastn or classifier.jar).
    @return: [str] The software absolute path.
    """
    path_locations = os.getenv("PATH").split(os.pathsep) + sys.path
    exec_path = None
    for current_location in path_locations:
        if exec_path is None and os.path.isfile(os.path.join(current_location, exec_name)):
            exec_path = os.path.abspath( os.path.join(current_location, exec_name) )
    if exec_path is None:
        raise Exception( "The software '" + exec_name + "' cannot be retrieved in path." )
    return exec_path


def prevent_shell_injections(argparse_namespace, excluded_args=None):
    """
    @summary: Raises an exception if one parameter contains a backquote or a semi-colon.
    @param argparse_namespace: [Namespase] The result of parser.parse_args().
    @param excluded_args: [list] List of unchecked parameters.
    """
    exceptions = list() if excluded_args is None else excluded_args
    for param_name in argparse_namespace.__dict__.keys():
        if not param_name in exceptions:
            param_val = getattr(argparse_namespace, param_name)
            if issubclass(param_val.__class__, list):
                new_param_val = list()
                for val in param_val:
                    if ';' in val.encode('utf8') or '`' in val.encode('utf8') or '|' in val.encode('utf8'):
                        raise Exception( "';' and '`' are unauthorized characters." ) 
            elif param_val is not None and issubclass(param_val.__class__, str):
                if ';' in param_val.encode('utf8') or '`' in param_val.encode('utf8') or '|' in param_val.encode('utf8'):
                    raise Exception( "';' and '`' are unauthorized characters." )


##################################################################################################################################################
#
# SEQIO
#
##################################################################################################################################################

class Sequence:
    def __init__(self, id, string, description=None, quality=None):
        """
        @param id : [str] Id of the sequence.
        @param string : [str] Sequence of the sequence.
        @param description : [str] The sequence description.
        @param quality : [str] The quality of the sequence (same length as string).
        """
        self.id = id
        self.description = description
        self.string = string
        self.quality = quality


class SequenceFileReader(object):
    @staticmethod
    def factory(filepath):
        if FastqIO.is_valid(filepath):
            return FastqIO(filepath)
        elif FastaIO.is_valid(filepath):
            return FastaIO(filepath)
        else:
            raise IOError( "The file " + filepath + " does not have a valid format for 'SequenceFileReader'." )


class FastqIO:
    def __init__(self, filepath, mode="r"):
        """
        @param filepath : [str] The filepath.
        @param mode : [str] Mode to open the file ('r', 'w', 'a').
        """
        self.filepath = filepath
        self.mode = mode
        if (mode in ["w", "a"] and filepath.endswith('.gz')) or (mode not in ["w", "a"] and is_gzip(filepath)):
            self.file_handle = gzip.open( filepath, mode )
        else:
            self.file_handle = open( filepath, mode )
        self.current_line_nb = 1
        self.current_line = None

    def __del__(self):
        self.close()

    def close(self):
        if hasattr(self, 'file_handle') and self.file_handle is not None:
            self.file_handle.close()
            self.file_handle = None
            self.current_line_nb = None

    def __iter__(self):
        seq_id = None
        seq_desc = None
        seq_str = None
        seq_qual = None
        try:
            for line in self.file_handle:
                line = line.rstrip()
                if (self.current_line_nb % 4) == 1:
                    fields = line[1:].split(None, 1)
                    seq_id = fields[0]
                    seq_desc = fields[1] if len(fields) == 2 else None
                elif (self.current_line_nb % 4) == 2:
                    seq_str = line
                elif (self.current_line_nb % 4) == 0:
                    seq_qual = line
                    yield Sequence( seq_id, seq_str, seq_desc, seq_qual )
                self.current_line_nb += 1
        except:
            raise IOError( "The line " + str(self.current_line_nb) + " in '" + self.filepath + "' cannot be parsed by " + self.__class__.__name__ + "." )

    def next_seq(self):
        """
        @summary : Returns the next sequence.
        @return : [Sequence] The next sequence.
        """
        seq_record = None
        try:
            # Header
            header = self.file_handle.readline().strip()
            fields = header[1:].split(None, 1)
            seq_id = fields[0]
            seq_desc = fields[1] if len(fields) == 2 else None
            self.current_line_nb += 1
            # Sequence
            seq_str = self.file_handle.readline().strip()
            self.current_line_nb += 1
            # Separator
            separator = self.file_handle.readline()
            self.current_line_nb += 1
            # Quality
            seq_qual = self.file_handle.readline().strip()
            self.current_line_nb += 1
            # Record
            seq_record = Sequence( seq_id, seq_str, seq_desc, seq_qual )
        except:
            raise IOError( "The line " + str(self.current_line_nb) + " in '" + self.filepath + "' cannot be parsed by " + self.__class__.__name__ + "." )
        return seq_record

    @staticmethod
    def is_valid(filepath):
        is_valid = False
        FH_in = FastqIO(filepath)
        try:
            seq_idx = 0
            previous = None
            while seq_idx < 10 and (seq_idx != 0 and previous is not None):
                previous = FH_in.next_seq()
                seq_idx += 1
            FH_in.close()
            # Cheack first header
            FH_in = FastqIO(filepath)
            if seq_idx == 0 or FH_in.file_handle.readline().startswith("@"):
                is_valid = True
        except:
            pass
        finally:
            FH_in.close()
        return is_valid

    def write(self, sequence_record):
        self.file_handle.write( self.seqToFastqLine(sequence_record) + "\n" )

    def seqToFastqLine(self, sequence):
        """
        @summary : Returns the sequence in fastq format.
        @param sequence : [Sequence] The sequence to process.
        @return : [str] The sequence.
        """
        seq = "@" + sequence.id + (" " + sequence.description if sequence.description is not None else "")
        seq += "\n" + sequence.string
        seq += "\n+"
        seq += "\n" + sequence.quality
        return seq


class FastaIO:
    def __init__(self, filepath, mode="r"):
        """
        @param filepath : [str] The filepath.
        @param mode : [str] Mode to open the file ('r', 'w', 'a').
        """
        self.filepath = filepath
        self.mode = mode
        if (mode in ["w", "a"] and filepath.endswith('.gz')) or (mode not in ["w", "a"] and is_gzip(filepath)):
            self.file_handle = gzip.open( filepath, mode )
        else:
            self.file_handle = open( filepath, mode )
        self.current_line_nb = 1
        self.current_line = None

    def __del__(self):
        self.close()

    def close(self):
        if hasattr(self, 'file_handle') and self.file_handle is not None:
            self.file_handle.close()
            self.file_handle = None
            self.current_line_nb = None

    def __iter__(self):
        seq_id = None
        seq_desc = None
        seq_str = None
        try:
            for line in self.file_handle:
                line = line.rstrip()
                self.current_line_nb += 1
                if line.startswith('>'):
                    if seq_id is not None:
                        seq_record = Sequence( seq_id, seq_str, seq_desc )
                        yield seq_record
                    # New seq
                    fields = line[1:].split(None, 1)
                    seq_id = fields[0]
                    seq_desc = fields[1] if len(fields) == 2 else None
                    seq_str = ""
                else:
                    seq_str += line
            if seq_id is not None:
                seq_record = Sequence( seq_id, seq_str, seq_desc )
                yield seq_record
        except:
            raise IOError( "The line " + str(self.current_line_nb) + " in '" + self.filepath + "' cannot be parsed by " + self.__class__.__name__ + "." )

    def next_seq(self):
        """
        @summary : Returns the next sequence.
        @return : [Sequence] The next sequence.
        """
        seq_record = None
        line = ""
        try:
            # First line in file
            if self.current_line_nb == 1:
                self.next_id = self.file_handle.readline().strip()
                self.current_line_nb += 1
            # Sequence
            seq_str = ""
            while not line.startswith('>'):
                seq_str += line.strip()
                line = self.file_handle.readline()
                if not line:
                    break
                self.current_line_nb += 1
            fields = self.next_id[1:].split(None, 1)
            seq_id = fields[0]
            seq_desc = fields[1] if len(fields) == 2 else None
            seq_record = Sequence( seq_id, seq_str, seq_desc )
            self.next_id = line # next seq_id
        except:
            raise IOError( "The line " + str(self.current_line_nb) + " in '" + self.filepath + "' cannot be parsed by " + self.__class__.__name__ + ".\n"
                            + "content : " + line )
        return seq_record

    @staticmethod
    def is_valid(filepath):
        is_valid = False
        FH_in = FastaIO(filepath)
        try:
            seq_idx = 0
            previous = None
            while seq_idx < 10 and (seq_idx != 0 and previous is not None):
                previous = FH_in.next_seq()
                seq_idx += 1
            is_valid = True
        except:
            pass
        finally:
            FH_in.close()
        return is_valid

    def write(self, sequence_record):
        self.file_handle.write( self.seqToFastaLine(sequence_record) + "\n" )

    def seqToFastaLine(self, sequence):
        """
        @summary : Returns the sequence in fasta format.
        @param sequence : [Sequence] The sequence to process.
        @return : [str] The sequence.
        """
        header = ">" + sequence.id + (" " + sequence.description if sequence.description is not None else "")
        return header + "\n" + sequence.string

##################################################################################################################################################
#
# CMD
#
##################################################################################################################################################

class Cmd:
    """
    @summary : Command wrapper.
    """
    def __init__(self, program, description, exec_parameters, version_parameters=None):
        """
        @param exec_parameters: [str] The parameters to execute the program. Two possibles syntaxes.
                                If the parameter contains the string '##PROGRAM##', this tag will be replaced by the program parameter before submit.
                                Otherwise the parameters will be added after the program in command line.
        @param version_parameters: [str] The parameters to get the program version. Two possibles syntaxes.
                                   If the parameter contains the string '##PROGRAM##', this tag will be replaced by the program parameter before submit.
                                   Otherwise the parameters will be added after the program in command line.
        """
        self.program = program
        self.description = description
        self.exec_parameters = exec_parameters
        self.version_parameters = version_parameters

    def get_cmd(self):
        """
        @summary : Returns the command line.
        @return : [str] The command line.
        """
        cmd = None
        if '##PROGRAM##' in self.exec_parameters:
            cmd = self.exec_parameters.replace('##PROGRAM##', self.program)
        else:
            cmd = self.program + ' ' + self.exec_parameters
        return cmd

    def get_version(self, location='stderr'):
        """
        @summary : Returns the program version number.
        @param location : [str] If the version command returns the version number on 'stdout' or on 'stderr'.
        @return : [str] version number if this is possible, otherwise this method return 'unknown'.
        """
        if self.version_parameters is None:
            return "unknown"
        else:
            try:
                cmd = self.program + ' ' + self.version_parameters
                if '##PROGRAM##' in self.exec_parameters:
                    cmd = self.version_parameters.replace('##PROGRAM##', self.program)
                p = Popen(cmd, shell=True, stdout=PIPE, stderr=PIPE)
                stdout, stderr = p.communicate()
                if location == 'stderr':
                    return stderr.strip()
                else:
                    return stdout.strip()
            except:
                raise Exception( "Version cannot be retrieve for the software '" + self.program + "'." )

    def parser(self, log_file):
        """
        @summary : Parse the command results to add information in log_file.
        @log_file : [str] Path to the sample process log file.
        """
        pass

    def submit(self, log_file=None):
        """
        @summary : Launch command, trace this action in log and parse results.
        @log_file : [str] Path to the sample process log file.
        """
        # Log
        if log_file is not None:
            FH_log = Logger( log_file )
            FH_log.write( '# ' + self.description + ' (' + os.path.basename(self.program) + ' version : ' + self.get_version() + ')\n' )
            FH_log.write( 'Command:\n\t' + self.get_cmd() + '\n' )
            FH_log.write( 'Execution:\n\tstart: ' + time.strftime("%d %b %Y %H:%M:%S", time.localtime()) + '\n' )
            FH_log.close()
        # Process
        subprocess.check_output( self.get_cmd(), shell=True )
        # Log
        if log_file is not None:
            FH_log = Logger( log_file )
            FH_log.write( '\tend:   ' + time.strftime("%d %b %Y %H:%M:%S", time.localtime()) + '\n' )
            FH_log.close()
            # Post-process results
            self.parser(log_file)

##################################################################################################################################################
#
# LOG AND TMP FILES
#
##################################################################################################################################################


class Logger:
    """
    @summary: Log file handler.
    """

    def __init__(self, filepath=None):
        """
        @param filepath: [str] The log filepath. [default : STDOUT]
        """
        self.filepath = filepath
        if self.filepath is not None and self.filepath is not sys.stdout:
            self.file_handle = open( self.filepath, "a" )
        else:
            self.file_handle = sys.stdout

    def __del__(self):
        """
        @summary: Closed file handler when the logger is detroyed.
        """
        self.close()

    def close(self):
        """
        @summary: Closed file handler.
        """
        if self.filepath is not None and self.filepath is not sys.stdout:
            if self.file_handle is not None:
                self.file_handle.close()
                self.file_handle = None

    def write(self, msg):
        """
        @summary: Writes msg on file.
        @param msg: [str] The message to write.
        """
        self.file_handle.write( msg )

    @staticmethod
    def static_write(filepath, msg):
        """
        @summary: Writes msg on file.
        @param filepath: [str] The log filepath. [default : STDOUT]
        @param msg: [str] The message to write.
        """
        if filepath is not None and filepath is not sys.stdout:
            FH_log = open( filepath, "a" )
            FH_log.write( msg )
            FH_log.close()
        else:
            sys.stdout.write( msg )


class TmpFiles:
    """
    @summary: Manager for temporary files.
    @note:
        tmpFiles = TmpFiles(out_dir)
        try:
            ...
            tmp_seq = tmpFiles.add( "toto.fasta" )
            ...
            tmp_log = tmpFiles.add( "log.txt" )
            ...
        finaly:
            tmpFiles.deleteAll()
    """
    def __init__(self, tmp_dir, prefix=None):
        """
        @param tmp_dir: [str] The temporary directory path.
        @param prefix: [str] The prefix added to each temporary file [default: <TIMESTAMP>_<PID>].
        """
        if prefix is None:
            prefix = str(time.time()) + "_" + str(os.getpid())
        self.files = list()
        self.dirs = list()
        self.tmp_dir = tmp_dir
        self.prefix = prefix

    def add(self, filename, prefix=None, dir=None):
        """
        @summary: Add a temporary file.
        @param filename: The filename without prefix.
        @param prefix: The prefix added [default: TmpFiles.prefix].
        @param dir: The directory path [default: TmpFiles.tmp_dir].
        @return: [str] The filepath.
        """
        # Default
        if prefix is None:
            prefix = self.prefix
        if dir is None:
            dir = self.tmp_dir
        # Process
        filepath = os.path.join(dir, prefix + "_" + filename)
        self.files.append(filepath)
        return filepath

    def add_dir(self, dirname, prefix=None, dir=None):
        """
        @summary: Add a temporary dir.
        @param filename: The dirname without prefix.
        @param prefix: The prefix added [default: TmpFiles.prefix].
        @param dir: The directory path [default: TmpFiles.tmp_dir].
        @return: [str] The filepath.
        """
        # Default
        if prefix is None:
            prefix = self.prefix
        if dir is None:
            dir = self.tmp_dir
        # Process
        dirpath = os.path.join(dir, prefix + "_" + dirname)
        self.dirs.append(dirpath)
        return dirpath

    def delete(self, filepath):
        """
        @summary: Deletes the specified temporary file.
        @param filepath: [str] The file path to delete.
        """
        self.files.remove(filepath)
        if os.path.exists(filepath): os.remove(filepath)

    def delete_dir(self, dirpath):
        """
        @summary: Deletes the specified temporary dir.
        @param filepath: [str] The file path to delete.
        """
        if dirpath in self.dirs: self.dirs.remove(dirpath)
            
        if os.path.exists(dirpath):
            for root, dirnames,filenames in os.walk(dirpath):
                for f in filenames:
                    if f in self.files: self.files.remove(os.path.join(dirpath,f))
                    if os.path.exists(os.path.join(dirpath,f)): os.remove(os.path.join(dirpath,f))
                for d in dirnames:
                    if d in self.dirs: self.dirs.remove(os.path.join(dirpath,d))
                    if os.path.exists(os.path.join(dirpath,d)): self.delete_dir(os.path.join(dirpath,d))
            os.rmdir(dirpath)

    def deleteAll(self):
        """
        @summary: Deletes all temporary files.
        """
        all_tmp_files = [tmp_file for tmp_file in self.files]
        for tmp_file in all_tmp_files:
            self.delete(tmp_file)
        
        all_tmp_dirs=[tmp_dir for tmp_dir in self.dirs]
        for tmp_dir in all_tmp_dirs:
            self.delete_dir(tmp_dir)
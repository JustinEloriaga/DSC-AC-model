function build_info = build_dsc_mex()
%BUILD_DSC_MEX Build the current constant-mean DSC native kernels.
% Uses the configured C++ compiler and MATLAB's bundled LAPACK/BLAS.
% On macOS, invoke Apple Command Line Tools directly, without mex creating
% temporary build directories. Every artifact uses the existing toolbox.
% The classical correlation kernel is built for regression/timing comparisons;
% the historical dynamic-mean FFBS kernel is not part of the current build.
% This helper does not change sampler defaults or compiler preferences.
% After building, add toolbox to the MATLAB path to call
% [ll,P,L,logdet,iterations,x] = dsc_corr_batch_mex_new(r,Z,mask,threads,guess)
% mu = dsc_constant_mean_mex(Y,h,P,mask,mu0,S0,z)
% threads is optional (default 1); native date workers are limited to 64.
root = fileparts(fileparts(mfilename('fullpath')));
target_dir = fullfile(root,'toolbox');
assert(isfolder(target_dir),'dsc:MexBuild','The toolbox directory must already exist.');
source = fullfile(target_dir,'dsc_corr_batch_mex.cpp');
new_source = fullfile(target_dir,'dsc_corr_batch_mex_new.cpp');
mean_source = fullfile(target_dir,'dsc_constant_mean_mex.cpp');
compiler = mex.getCompilerConfigurations('C++','Selected');
use_command_line_tools = ismac;
if isempty(compiler) && ~use_command_line_tools
    error('dsc:MexCompiler', ...
        'No C++ MEX compiler is configured. Run mex -setup C++ and retry build_dsc_mex.');
end
fprintf('Constant-mean DSC MEX build: MATLAB %s (%s), %s\n', ...
    version,version('-release'),computer('arch'));
fprintf('Source: %s\nOutput: %s\n',source,target_dir);
fprintf('Libraries: MATLAB bundled mwlapack and mwblas; optional std::thread workers, no OpenMP.\n');
started = tic;
clear dsc_corr_batch_mex dsc_corr_batch_mex_new dsc_constant_mean_mex
artifact = fullfile(target_dir,['dsc_corr_batch_mex.' mexext]);
new_artifact = fullfile(target_dir,['dsc_corr_batch_mex_new.' mexext]);
mean_artifact = fullfile(target_dir,['dsc_constant_mean_mex.' mexext]);
if use_command_line_tools
    [compiler_name,compiler_version] = build_with_command_line_tools(source,artifact);
    build_with_command_line_tools(new_source,new_artifact);
    build_with_command_line_tools(mean_source,mean_artifact);
    build_method = 'Apple Command Line Tools (direct compiler invocation)';
else
    compiler = compiler(1);
    compiler_name = compiler.Name; compiler_version = compiler.Version;
    fprintf('Compiler: %s, version %s\n',compiler_name,compiler_version);
    threading_flags = {};
    if isunix && ~ismac
        threading_flags = {'CXXFLAGS=$CXXFLAGS -pthread','LDFLAGS=$LDFLAGS -pthread'};
    end
    mex('-R2018a','-O','-v',threading_flags{:},source,'-lmwlapack','-lmwblas', ...
        '-outdir',target_dir,'-output','dsc_corr_batch_mex');
    mex('-R2018a','-O','-v',threading_flags{:},new_source,'-lmwlapack','-lmwblas', ...
        '-outdir',target_dir,'-output','dsc_corr_batch_mex_new');
    mex('-R2018a','-O','-v',threading_flags{:},mean_source,'-lmwlapack','-lmwblas', ...
        '-outdir',target_dir,'-output','dsc_constant_mean_mex');
    build_method = 'MATLAB mex';
end
build_info = struct('artifact',artifact,'source',source, ...
    'new_artifact',new_artifact,'new_source',new_source, ...
    'mean_artifact',mean_artifact,'mean_source',mean_source, ...
    'matlab_version',version,'matlab_release',version('-release'), ...
    'architecture',computer('arch'),'compiler_name',compiler_name, ...
    'compiler_version',compiler_version,'build_method',build_method,'seconds',toc(started), ...
    'libraries',{{'mwlapack','mwblas'}},'openmp',false, ...
    'native_threads_default',1,'native_threads_limit',64, ...
    'model','constant_mu_v1','mean_backend','constant_gaussian', ...
    'legacy_correlation_reference',artifact,'dynamic_mean_kernel_built',false);
fprintf('Built %s, %s and %s in %.2f seconds.\n',artifact,new_artifact,mean_artifact,build_info.seconds);
end

function [name,compiler_version] = build_with_command_line_tools(source,artifact)
[status,compiler_version] = system('/usr/bin/xcrun clang++ --version');
if status ~= 0
    error('dsc:MexCompiler','Neither a selected MEX compiler nor Apple Command Line Tools is available.');
end
name = 'Apple Command Line Tools clang++';
compiler_version = strtrim(compiler_version);
fprintf('Compiler: %s\n',compiler_version);
switch computer('arch')
    case 'maca64', architecture = 'arm64'; deployment = '-mmacosx-version-min=12.0';
    case 'maci64', architecture = 'x86_64'; deployment = '-mmacosx-version-min=10.15';
    otherwise, error('dsc:MexCompiler','Unsupported macOS MATLAB architecture.');
end
% Match the interleaved-complex C Matrix API selected by mex -R2018a.
% MATLAB supplies the API-version shim and symbol map as well as libraries.
library_dir = fullfile(matlabroot,'bin',computer('arch'));
export_map = fullfile(matlabroot,'extern','lib',computer('arch'), ...
    'c_exportsmexfileversion.map');
arguments = {'/usr/bin/xcrun','clang++','-std=c++14','-O2','-DNDEBUG', ...
    '-fno-common','-fwrapv','-ffp-contract=off','-fexceptions', ...
    '-arch',architecture,deployment,'-bundle','-stdlib=libc++', ...
    '-DMATLAB_MEX_FILE','-DMATLAB_DEFAULT_RELEASE=R2018a', ...
    ['-I' fullfile(matlabroot,'extern','include')],source, ...
    fullfile(matlabroot,'extern','version','cpp_mexapi_version.cpp'), ...
    ['-L' library_dir],'-lmx','-lmex','-lmat','-lmwlapack','-lmwblas', ...
    ['-Wl,-rpath,' library_dir],['-Wl,-exported_symbols_list,' export_map], ...
    '-o',artifact};
command = strjoin(cellfun(@shell_quote,arguments,'UniformOutput',false),' ');
[status,message] = system(command,'-echo');
if status ~= 0, error('dsc:MexBuild','Native correlation build failed: %s',message); end
end

function quoted = shell_quote(value)
apostrophe = char(39);
quoted = [apostrophe strrep(char(value),apostrophe, ...
    [apostrophe '"' apostrophe '"' apostrophe]) apostrophe];
end

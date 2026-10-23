function build_info = build_constant_mean_mex()
%BUILD_CONSTANT_MEAN_MEX Build the exact constant Gaussian mean kernel.
% Uses MATLAB's bundled LAPACK/BLAS and an existing toolbox directory.
root=fileparts(fileparts(mfilename('fullpath')));
target=fullfile(root,'toolbox');
assert(isfolder(target),'dsc:MexBuild','The existing toolbox directory is required.');
source=fullfile(target,'dsc_constant_mean_mex.cpp');
artifact=fullfile(target,['dsc_constant_mean_mex.' mexext]);
compiler=mex.getCompilerConfigurations('C++','Selected');
clear dsc_constant_mean_mex
started=tic;
if ismac
    [compiler_name,compiler_version]=build_with_command_line_tools(source,artifact);
    method='Apple Command Line Tools (direct compiler invocation)';
elseif isempty(compiler)
    error('dsc:MexCompiler','No C++ MEX compiler is configured. Run mex -setup C++ and retry.');
else
    compiler=compiler(1);
    compiler_name=compiler.Name; compiler_version=compiler.Version;
    mex('-R2018a','-O',source,'-lmwlapack','-lmwblas', ...
        '-outdir',target,'-output','dsc_constant_mean_mex');
    method='MATLAB mex';
end
build_info=struct('artifact',artifact,'source',source,'matlab_version',version, ...
    'matlab_release',version('-release'),'architecture',computer('arch'), ...
    'compiler_name',compiler_name,'compiler_version',compiler_version, ...
    'build_method',method,'seconds',toc(started),'libraries',{{'mwlapack','mwblas'}});
fprintf('Built %s in %.2f seconds.\n',artifact,build_info.seconds);
end

function [name,compiler_version]=build_with_command_line_tools(source,artifact)
[status,compiler_version]=system('/usr/bin/xcrun clang++ --version');
if status~=0
    error('dsc:MexCompiler','Neither a selected MEX compiler nor Apple Command Line Tools is available.');
end
name='Apple Command Line Tools clang++'; compiler_version=strtrim(compiler_version);
switch computer('arch')
    case 'maca64', architecture='arm64'; deployment='-mmacosx-version-min=12.0';
    case 'maci64', architecture='x86_64'; deployment='-mmacosx-version-min=10.15';
    otherwise, error('dsc:MexCompiler','Unsupported macOS MATLAB architecture.');
end
library_dir=fullfile(matlabroot,'bin',computer('arch'));
export_map=fullfile(matlabroot,'extern','lib',computer('arch'),'c_exportsmexfileversion.map');
arguments={'/usr/bin/xcrun','clang++','-std=c++14','-O2','-DNDEBUG', ...
    '-fno-common','-fwrapv','-ffp-contract=off','-fexceptions', ...
    '-arch',architecture,deployment,'-bundle','-stdlib=libc++', ...
    '-DMATLAB_MEX_FILE','-DMATLAB_DEFAULT_RELEASE=R2018a', ...
    ['-I' fullfile(matlabroot,'extern','include')],source, ...
    fullfile(matlabroot,'extern','version','cpp_mexapi_version.cpp'), ...
    ['-L' library_dir],'-lmx','-lmex','-lmat','-lmwlapack','-lmwblas', ...
    ['-Wl,-rpath,' library_dir],['-Wl,-exported_symbols_list,' export_map],'-o',artifact};
command=strjoin(cellfun(@shell_quote,arguments,'UniformOutput',false),' ');
[status,message]=system(command,'-echo');
if status~=0, error('dsc:MexBuild','Native constant mean build failed: %s',message); end
end

function quoted=shell_quote(value)
apostrophe=char(39);
quoted=[apostrophe strrep(char(value),apostrophe,[apostrophe '"' apostrophe '"' apostrophe]) apostrophe];
end

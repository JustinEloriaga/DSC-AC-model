function info=build_correlation_new_mex()
%BUILD_CORRELATION_NEW_MEX Build only the Newton correlation kernel in place.
% No directory creation, compiler preference changes, or sampler execution.
root=fileparts(fileparts(mfilename('fullpath')));
source=fullfile(root,'toolbox','dsc_corr_batch_mex_new.cpp');
artifact=fullfile(root,'toolbox',['dsc_corr_batch_mex_new.' mexext]);
compiler=mex.getCompilerConfigurations('C++','Selected');
started=tic;
clear dsc_corr_batch_mex_new
if ismac
    [status,compiler_version]=system('/usr/bin/xcrun clang++ --version');
    assert(status==0,'dsc:MexCompiler','Apple Command Line Tools are unavailable.');
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
    assert(status==0,'dsc:MexBuild','Correlation build failed: %s',message);
    compiler_name='Apple Command Line Tools clang++';
    build_method='Direct compiler invocation';
else
    assert(~isempty(compiler),'dsc:MexCompiler','No C++ MEX compiler is configured.');
    flags={};
    if isunix&&~ismac, flags={'CXXFLAGS=$CXXFLAGS -pthread','LDFLAGS=$LDFLAGS -pthread'}; end
    mex('-R2018a','-O',flags{:},source,'-lmwlapack','-lmwblas', ...
        '-outdir',fileparts(artifact),'-output','dsc_corr_batch_mex_new');
    compiler_name=compiler(1).Name; compiler_version=compiler(1).Version;
    build_method='MATLAB mex';
end
info=struct('source',source,'artifact',artifact,'seconds',toc(started), ...
    'compiler',compiler_name,'compiler_version',strtrim(compiler_version), ...
    'build_method',build_method,'matlab_version',version,'architecture',computer('arch'));
disp(info);
end

function quoted=shell_quote(value)
apostrophe=char(39);
quoted=[apostrophe strrep(char(value),apostrophe, ...
    [apostrophe '"' apostrophe '"' apostrophe]) apostrophe];
end

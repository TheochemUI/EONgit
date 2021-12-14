# Generates a development environment for eonclient
#
# Usage examples:
#
# To create a dev environment:
#
#   nix-shell
#
# With clang:
#
#   nix-shell --argstr doomdir clang
#
{ withEonclient ? false, compiler ? "gcc" }:
let
  sources = import ./nix/sources.nix;
  pkgs = import sources.nixpkgs {};
  mach-nix = import (
    builtins.fetchGit {
      url = "https://github.com/DavHau/mach-nix.git";
      ref = "refs/tags/3.3.0";
    }
  ) {
    pkgs = pkgs;
    python = "python38";
  };
  customPython = mach-nix.mkPython {
    requirements = builtins.readFile ./requirements.txt;
    packagesExtra = [
      # "https://github.com/psf/requests/tarball/2a7832b5b06d"   # from tarball url
      ./. # from local path
      # mach-nix.buildPythonPackage { ... };                     # from package
    ];
  };
  compilerEnv = (
    if compiler == "gcc" then pkgs.gcc10Stdenv
    else if compiler == "clang" then pkgs.clang10Stdenv
    else pkgs.stdenv
  );
  mkShellNewEnv = pkgs.mkShell.override { stdenv = compilerEnv; };
  eigen339 = pkgs.eigen.overrideAttrs (
    old: rec {
      version = "3.3.9";
      stdenv = compilerEnv;
      src = pkgs.fetchFromGitLab {
        owner = "libeigen";
        repo = "eigen";
        rev = "${version}";
        sha256 = "0m4h9fd5s1pzpncy17r3w0b5a6ywqjajmnr720ndb7fc4bn0dhi4";
      };
      # From https://github.com/foolnotion/aoc2020/blob/master/eigen_include_dir.patch
      patches = [ ./eigen_include_dir.patch ];
    }
  );
  macHook = ''
    # eonclient
    export PATH=$(pwd)/client/builddir:$PATH
  '';
  linuxHook = ''
    # eonclient
    export PATH=$(pwd)/client/builddir:$PATH
    # Locale
    export LOCALE_ARCHIVE=${pkgs.glibcLocales}/lib/locale/locale-archive
  '';
  myCmop = (
    if compiler == "gcc" then (
      pkgs.wrapCC (
        pkgs.gcc10.cc.override {
          langFortran = true;
          langCC = true;
          langC = true;
          enableShared = true;
          enableMultilib = false;
          staticCompiler = false;
          profiledCompiler = false;
          hardeningEnable = [ "pic" ];
        }
      )
    )
    else if compiler == "clang" then (
      compilerEnv
    ) else pkgs.stdEnv
  );
  eonclient = pkgs.callPackage ./default.nix {};
in
mkShellNewEnv {
  nativeBuildInputs = with pkgs; [ cmake blas (if compiler == "gcc" then myCmop myCmop.cc.lib else myCmop) openblas ninja ];
  buildInputs = with pkgs; [
    gtest
    bashInteractive
    which
    customPython
    ninja
    meson
    myCmop
    #  gcc10Stdenv
    #  gfortran
    #   valgrind
    (if pkgs.stdenv.isDarwin then null else lldb)
    (if withEonclient then (if pkgs.stdenv.isDarwin then null else eonclient) else null)
    (if compiler == "clang" then gfortran else null)
    graphviz
    # gfortran
    # gfortran.cc
    # myCmop.cc

    zstd
    zlib
    lzma
    bzip2
    openblas

    fmt
    abseil-cpp
    boost175
    eigen339
  ];
  shellHook = if pkgs.stdenv.isDarwin then macHook else linuxHook;
}

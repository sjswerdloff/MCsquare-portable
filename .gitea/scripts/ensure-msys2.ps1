# Ensure a verified MSYS2 with UCRT64 gcc + libgomp on this Windows runner, at a VERSIONED path
# (C:\mcsq-win\msys64-<release>). Idempotent; called by ci-windows.yml and windows-vs-linux-field-edge.yml,
# because the `windows` label can land a job on either box. Throws on any failure.
#
# NON-DESTRUCTIVE BY DESIGN: this script never deletes or modifies an existing tree it did not create
# at this path. Older trees (e.g. C:\mcsq-win\msys64 from earlier, unverified runs) are left alone;
# removing them is sjswerdloff's decision. If this path exists without the marker proving it was
# unpacked from the pinned installer, the script fails closed and changes nothing.
$ErrorActionPreference = 'Continue'  # native stderr is not an error; results are checked explicitly
"$env:COMPUTERNAME  $((Get-CimInstance Win32_Processor).Name)"
$ver  = '20260927'
$want = 'ad336cccfda47758b5e15cda993fbba421115cb0b126697daef1ee4dfe37209f'  # GitHub release digest = repo.msys2.org file
$m    = "C:\mcsq-win\msys64-$ver"
$marker = "$m\.installed-from-sha256"
New-Item -ItemType Directory -Force C:\mcsq-win | Out-Null

if (Test-Path $m) {
  $installed = if (Test-Path $marker) { (Get-Content $marker -Raw).Trim() } else { '' }
  if ($installed -ne $want) { throw "$m exists but was not unpacked from the pinned installer (marker '$installed'): not touching it; ask sjswerdloff" }
} else {
  $stage = "C:\mcsq-win\msys2-stage-$ver-$env:GITHUB_RUN_ID-$PID"   # unique: a cancelled run's leftover never blocks
  $sfx = "C:\mcsq-win\msys2-base-x86_64-$ver.sfx.exe"
  if (-not (Test-Path $sfx)) {
    curl.exe -fsSL -o $sfx "https://repo.msys2.org/distrib/x86_64/msys2-base-x86_64-$ver.sfx.exe"
    if ($LASTEXITCODE) { throw "MSYS2 download failed" }
  }
  $got = (Get-FileHash $sfx -Algorithm SHA256).Hash.ToLower()
  if ($got -ne $want) { throw "MSYS2 installer sha256 $got, expected ${want}: not run (file left at $sfx for inspection)" }
  # The archive unpacks to <dir>\msys64; unpack into a fresh staging dir and move that into place.
  & $sfx -y "-o$stage" | Out-Null
  if (-not (Test-Path "$stage\msys64\usr\bin\bash.exe")) { throw "MSYS2 did not unpack" }
  # Marker first, then one rename: the tree and its marker appear at $m together, so a cancelled job
  # (a push cancels in-flight runs) can leave only a staging dir, never an unmarked tree at $m.
  Set-Content -NoNewline "$stage\msys64\.installed-from-sha256" $want
  Move-Item "$stage\msys64" $m
  Remove-Item $stage   # the now-empty staging dir this script created
}

$env:MSYSTEM = 'UCRT64'; $env:CHERE_INVOKING = '1'
if (-not (Test-Path "$m\ucrt64\bin\gcc.exe") -or -not (Test-Path "$m\ucrt64\lib\libgomp.dll.a")) {
  # First run initialises the keyring; a core update can end the shell, so update twice.
  & "$m\usr\bin\bash.exe" -lc 'pacman -Syuu --noconfirm' 2>&1 | Select-Object -Last 3
  & "$m\usr\bin\bash.exe" -lc 'pacman -Syuu --noconfirm' 2>&1 | Select-Object -Last 3
  & "$m\usr\bin\bash.exe" -lc 'pacman -S --needed --noconfirm mingw-w64-ucrt-x86_64-gcc mingw-w64-ucrt-x86_64-libgomp make' 2>&1 | Select-Object -Last 3
  if (-not (Test-Path "$m\ucrt64\bin\gcc.exe") -or -not (Test-Path "$m\ucrt64\lib\libgomp.dll.a")) { throw "UCRT64 gcc/libgomp not installed" }
}
"MSYS2 ready at $m"

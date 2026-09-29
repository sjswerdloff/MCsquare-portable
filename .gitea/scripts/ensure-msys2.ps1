# Ensure a verified MSYS2 (pinned dated release, SHA-256 checked) with UCRT64 gcc + libgomp at
# C:\mcsq-win\msys64 on this Windows runner. Idempotent; called by ci-windows.yml and
# windows-vs-linux-field-edge.yml, because the `windows` label can land a job on either box.
# Throws on any failure, so the calling step fails.

$ErrorActionPreference = 'Continue'  # native stderr is not an error; results are checked explicitly
"$env:COMPUTERNAME  $((Get-CimInstance Win32_Processor).Name)"
$m = 'C:\mcsq-win\msys64'
$ver = '20260927'
$want = 'ad336cccfda47758b5e15cda993fbba421115cb0b126697daef1ee4dfe37209f'
$sfx = "C:\mcsq-win\msys2-base-x86_64-$ver.sfx.exe"
New-Item -ItemType Directory -Force C:\mcsq-win | Out-Null
# Earlier runs unpacked an unverified "latest" installer. Keep that tree only if its installer
# hashes to the pinned digest; otherwise remove it and start again from a verified download.
$marker = "$m\.installed-from-sha256"
$installed = if (Test-Path $marker) { (Get-Content $marker -Raw).Trim() } else { '' }
if ($installed -ne $want -and (Test-Path "$m\usr\bin\bash.exe")) {
  $old = 'C:\mcsq-win\msys2.sfx.exe'
  $oldsha = if (Test-Path $old) { (Get-FileHash $old -Algorithm SHA256).Hash.ToLower() } else { 'none' }
  if ($oldsha -eq $want) { "existing tree came from the pinned installer: kept"; Set-Content -NoNewline $marker $want; $installed = $want }
  else { "existing tree came from installer sha256 $oldsha, not the pinned one: removing it"; Remove-Item -Recurse -Force $m }
}
if (-not (Test-Path "$m\usr\bin\bash.exe")) {
  curl.exe -fsSL -o $sfx "https://repo.msys2.org/distrib/x86_64/msys2-base-x86_64-$ver.sfx.exe"
  if ($LASTEXITCODE) { throw "MSYS2 download failed" }
  $got = (Get-FileHash $sfx -Algorithm SHA256).Hash.ToLower()
  if ($got -ne $want) { Remove-Item -Force $sfx; throw "MSYS2 installer sha256 $got, expected ${want}: not run" }
  & $sfx -y -oC:\mcsq-win\ | Out-Null
  if (-not (Test-Path "$m\usr\bin\bash.exe")) { throw "MSYS2 did not unpack" }
  Set-Content -NoNewline $marker $want
}
# Undo any safe.directory=* an earlier version of this workflow may have set in MSYS2's git.
& "$m\usr\bin\bash.exe" -lc 'command -v git >/dev/null && git config --global --unset-all safe.directory "^[*]$" ; true' 2>&1 | Out-Null
if (-not (Test-Path "$m\ucrt64\bin\gcc.exe")) {
  $env:MSYSTEM = 'UCRT64'; $env:CHERE_INVOKING = '1'
  # First run initialises the keyring; a core update can end the shell, so update twice.
  & "$m\usr\bin\bash.exe" -lc 'pacman -Syuu --noconfirm' 2>&1 | Select-Object -Last 5
  & "$m\usr\bin\bash.exe" -lc 'pacman -Syuu --noconfirm' 2>&1 | Select-Object -Last 5
  & "$m\usr\bin\bash.exe" -lc 'pacman -S --needed --noconfirm mingw-w64-ucrt-x86_64-gcc mingw-w64-ucrt-x86_64-libgomp make' 2>&1 | Select-Object -Last 5
  if (-not (Test-Path "$m\ucrt64\bin\gcc.exe")) { throw "UCRT64 gcc not installed" }
}
# OpenMP is a separate package (GCC 16 in MSYS2); a box that got gcc before this line may lack it.
if (-not (Test-Path "$m\ucrt64\lib\libgomp.dll.a")) {
  $env:MSYSTEM = 'UCRT64'
  & "$m\usr\bin\bash.exe" -lc 'pacman -S --needed --noconfirm mingw-w64-ucrt-x86_64-libgomp' 2>&1 | Select-Object -Last 3
  if (-not (Test-Path "$m\ucrt64\lib\libgomp.dll.a")) { throw "libgomp not installed" }
}


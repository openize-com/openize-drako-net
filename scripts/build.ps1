$EnvScript = ".build.env.ps1"
$OutputDir = "..\output"




if(Test-Path $EnvScript)
{
    . $EnvScript
}
if("${Env:TIMESTAMPER_URL}" -eq "")
{
    Write-Error "ERROR: Env:TIMESTAMPER_URL must points to a valid timestamp service url"
    return
}
if("${Env:CERT_SUBJECT}" -eq "")
{
    Write-Error "ERROR: Env:CERT_SUBJECT must points to a valid certificate subject name"
    return
}

$Version=Select-Xml -Path ..\Openize.Drako\Openize.Drako.csproj -XPath '/Project/PropertyGroup/Version' | ForEach-Object { $_.Node.InnerXML }
$nupkg="bin\Release\Openize.Drako.$Version.nupkg"

Write-Host "Building Openize.Drako Version $Version"
try
{
    Push-Location ../Openize.Drako
    dotnet clean
    if(Test-Path $nupkg) {
        Remove-Item $nupkg
    }



    dotnet pack /p:SignTool=${Env:SIGN_TOOL} /p:CertFingerprint=${Env:CERT_FINGERPRINT} /p:TimestamperUrl=${Env:TIMESTAMPER_URL}
    if($LASTEXITCODE -eq '1') {
        Write-Host "Build failed."
        return
    }
    Write-Host "Signing the NuGet package"
    dotnet nuget sign $nupkg --certificate-subject-name "${Env:CERT_SUBJECT}" --timestamper ${Env:TIMESTAMPER_URL}/legacy

    if(!(Test-Path $OutputDir))
    {
        New-Item -Type Directory $OutputDir > $null
    }
    $Output="$OutputDir\$(Split-Path $nupkg -Leaf)"
    if(Test-Path $Output)
    {
        Remove-Item $Output
    }
    Move-Item $nupkg $Output
    Write-Host "Output file: $Output"

}
finally
{
    Pop-Location
}
$clean_dirs = @( "..\Openize.Drako\bin", "..\Openize.Drako\obj", "..\Openize.Drako.App\bin", "..\Openize.Drako.App\obj", "..\Openize.Drako.Tests\bin", "..\Openize.Drako.Tests\obj")
foreach($dir in $clean_dirs)
{
    if(Test-Path $dir)
    {
        Remove-Item -Recurse -Force $dir
    }
}


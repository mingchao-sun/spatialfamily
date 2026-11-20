# Local preview of the spatialfamily ribbon icon. Run in an ADMINISTRATOR PowerShell:
#   powershell -ExecutionPolicy Bypass -File dev\icon_preview.ps1            (apply)
#   powershell -ExecutionPolicy Bypass -File dev\icon_preview.ps1 -Restore   (remove only this icon's rule)
#   powershell -ExecutionPolicy Bypass -File dev\icon_preview.ps1 -RestoreAll (put the original stylesheet back)
# It adds one CSS rule to jamovi's client stylesheet; the file is backed up once as .orig.
# Each icon's rule is tagged /*icon:NAME*/ so that scripts of several modules can coexist:
# applying one icon never removes the rules of the others.
# The ribbon button's data-name is the lower-cased menu group name, here "spatial" (menuGroup: Spatial).
# Close all jamovi windows first, and start jamovi again afterwards. A jamovi update overwrites the change.
param([switch]$Restore, [switch]$RestoreAll)
$name = "spatial"
$assets = "C:\Program Files\jamovi 28.17.0.0\Resources\client\assets"
$css = Get-ChildItem $assets -Filter "main-*.css" | Select-Object -First 1
if (-not $css) { throw "main-*.css not found in $assets" }
$bak = "$($css.FullName).orig"
$utf8 = New-Object Text.UTF8Encoding($false)
if ($RestoreAll) {
    if (Test-Path $bak) { Copy-Item $bak $css.FullName -Force; Remove-Item $bak; "Restored $($css.Name)" } else { "No backup found, nothing to restore" }
    return
}
if (-not (Test-Path $bak)) { Copy-Item $css.FullName $bak }
$text = [IO.File]::ReadAllText($css.FullName)
# drop an earlier rule of this icon (and nothing else)
$text = [regex]::Replace($text, "\r?\n/\*icon:$name\*/[^\r\n]*", "")
if ($Restore) {
    [IO.File]::WriteAllText($css.FullName, $text, $utf8)
    "Removed the $name rule from $($css.Name)"
    return
}
$svgPath = Join-Path (Split-Path $PSScriptRoot -Parent) "icon\analysis-spatial.svg"
$b64 = [Convert]::ToBase64String([IO.File]::ReadAllBytes($svgPath))
$rule = "`n/*icon:$name*/.jmv-ribbon-button[data-name=$name] .jmv-ribbon-button-icon{background-image:url(`"data:image/svg+xml;base64,$b64`")}"
[IO.File]::WriteAllText($css.FullName, $text.TrimEnd() + $rule + "`n", $utf8)
"Rule '$name' added to $($css.Name). Restart jamovi to see the icon."

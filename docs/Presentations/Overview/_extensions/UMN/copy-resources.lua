-- copy resources hack until Quarto devs fix the issue

local function copyResource(file)

  path = quarto.utils.resolvePath(file)
  
  quarto.doc.addFormatResource(path)

end

function Header(el)

  copyResource('D2D-maroon-blk.png')
  copyResource('M_gold.png')

  return el

end

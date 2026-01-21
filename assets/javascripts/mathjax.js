window.MathJax = {
    tex: {
        inlineMath: [["\\(", "\\)"], ["$", "$"]],
        displayMath: [["\\[", "\\]"], ["$$", "$$"]],
        processEscapes: true,
        processEnvironments: true
    },
    options: {
        // Process math in the main content area
        // ignoreHtmlClass and processHtmlClass removed to allow processing everywhere
        renderActions: {
            // Fix HTML entities in LaTeX before processing
            fixHtmlEntities: [5, (doc) => {
                for (const math of doc.math) {
                    math.math = math.math
                        .replace(/&amp;/g, '&')
                        .replace(/&lt;/g, '<')
                        .replace(/&gt;/g, '>')
                        .replace(/&quot;/g, '"');
                }
            }, '']
        }
    },
    svg: {
        fontCache: "global"
    }
};

document$.subscribe(() => {
    MathJax.typesetClear()
    MathJax.texReset()
    MathJax.typesetPromise()
})
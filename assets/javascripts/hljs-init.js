document.addEventListener('DOMContentLoaded', function () {
    // Initialize highlight.js
    hljs.highlightAll();

    // Optional: Add language labels
    document.querySelectorAll('pre code').forEach((block) => {
        const language = block.className.match(/language-(\w+)/);
        if (language) {
            const label = document.createElement('div');
            label.className = 'hljs-language-label';
            label.textContent = language[1];
            block.parentNode.insertBefore(label, block);
        }
    });
});